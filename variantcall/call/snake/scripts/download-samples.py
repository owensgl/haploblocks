#!/usr/bin/env python3

"""
  Usage: download-samples [OPT]* OUTPUTDIR

  Reads a list of URLS from standard input, downloads them,
  verifies checksums (if available), and places them in OUTPUTDIR,
  using content digests as filenames.

  Files that already exist are skipped.

  Each input line must use the syntax:

     URL [md5 [MD5URL]]
     ...

"""

import argparse
import yaml
import logging
import sys
import os, os.path
import errno
import requests
import tempfile
import shutil
import urllib.parse
import hashlib
import time
import concurrent.futures

class DownloadError(Exception):
    pass

class ContentCache(object):
    def __init__(self, cachedir):
        self.cachedir = cachedir

    def _cache_entry(self, digest_type, digest):
        return os.path.join(self.cachedir, digest_type.upper() + "_" + digest.strip().lower())

    def get(self, digest_type, digest):
        if not digest_type or not digest:
            return None

        link_name = self._cache_entry(digest_type, digest)
        try:
            target = os.readlink(link_name)
            if os.path.exists(link_name):
                return os.path.realpath(link_name)
            return None
        except OSError as ose:
            if ose.errno == errno.ENOENT:
                return None
            raise

    def put(self, file_path, digest_type, digest):
        link_name = self._cache_entry(digest_type, digest)
        basename = os.path.basename(file_path)
        dst_name = os.path.join(self.cachedir, basename)

        if os.path.exists(dst_name):
            raise DownloadError("Naming collision: %s is already present. Consider renaming." % (dst_name,))

        with tempfile.NamedTemporaryFile(prefix="tmp.digest.", dir=self.cachedir, delete=False) as digest_fd:
            digest_fd.write((digest + " " + basename + "\n").encode('utf-8'))
            os.rename(digest_fd.name, dst_name + "." + digest_type.lower())

        os.rename(file_path, dst_name)
        try:
            os.symlink(basename, link_name)
        except OSError as ose:
            target = os.readlink(link_name)
            if target != basename:
                raise DownloadError("Digest %s => %s should point to %s" % (link_name, target, basename))

        return dst_name

def setup_logging():
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stderr)

    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    root.addHandler(ch)

log = logging.getLogger("download")

def load_creds(credsfile):
    with open(credsfile, "r") as stream:
        try:
            creds = yaml.load(stream)
            return creds
        except yaml.YAMLError as exc:
            raise

def _do_download(dl_dir, url, creds, digest_obj=None, logprefix=""):
    log.info("%sdownloading %s ==> %s (username=%s)", logprefix, url, dl_dir, creds['username'])
    if creds['username']:
        r = requests.post(url, data={'j_username': creds['username'],
                                     'j_password': creds['password']},
                          stream=True)
    else:
        r = requests.get(url, stream=True)

    if len(r.history) > 0:
        oldest = r.history[0]
        if oldest.status_code > 300 and oldest.status_code < 400:
            log.info("%s  redirected %s %s", logprefix, oldest.status_code, oldest.headers.get("location"))
            if "j_security_check" in oldest.headers.get("location", ""):
                log.info("%s  download failed. bad credentials.", logprefix)
                raise Exception("download failed")

    content_length = r.headers.get('content-length', None)
    if content_length:
        content_length = int(content_length, 10)

    log.info("%s  content-Length %s", logprefix, content_length)

    #for hdr in r.headers:
    #    log.debug("%s  %s: %s", logprefix, hdr, r.headers[hdr])

    raw_data = r.raw

    download_output = os.path.join(dl_dir, urllib.parse.unquote(os.path.basename(url)))

    siz = 0
    chksiz = 256*1024
    last_update = 0
    with open(download_output, "wb") as fd:
        chunk = 'x'
        while chunk:
            chunk = raw_data.read(chksiz)
            siz += len(chunk)
            if digest_obj:
                digest_obj.update(chunk)
            fd.write(chunk)
            if siz % (5*1024*1024) < chksiz:
                now = time.time()
                if last_update + 5.0 < now:
                    last_update = now
                    prog_siz = siz // (1024*1024)
                    if content_length:
                        prog_pct = siz / (content_length + 1) * 100
                        log.info("%s  progress %6d / %d MiB (%5.2f%%)", logprefix,
                                 prog_siz, (content_length + 512*1024) // (1024*1024), prog_pct)
                    else:
                        log.info("%s  progress %6d / -- MiB   --", logprefix, prog_siz)
    log.info("%s  done. file %s written.", logprefix, download_output)

    return download_output

def _extract_digest_entry(digest_file, entry_name):
    with open(digest_file, "r") as fd:
        for line in fd:
            line = line.strip()
            if not line:
                continue
            dig, nam = line.split(maxsplit=1)
            if nam == entry_name:
                return dig
    return None

def download_content(cache, work_dir, url, digest_type, digest_url, creds, serialno):
    tempdir = tempfile.mkdtemp(prefix=".download-samples-", dir=work_dir)
    basename = os.path.basename(url)
    expected_digest = None
    digest_obj = None

    try:
        if digest_type:
            digest_type = digest_type.upper()
        if digest_type in ('MD5',):
            md5_file = _do_download(tempdir, digest_url, creds,
                                    logprefix="%04d.md5 " % (serialno,))

            expected_digest = _extract_digest_entry(md5_file, basename)
            if not expected_digest:
                raise Exception("cannot obtain expected digest for entry: %s" % (basename,))

            digest_obj = hashlib.md5()

        elif digest_type:
            log.error("%04d digest type %s not supported", serialno, digest_type)
            raise Exception("cannot download content")

        log.info("%04d looking up content cache for %s:%s", serialno, digest_type, expected_digest)

        cached_path = cache.get(digest_type, expected_digest)
        if cached_path:
            log.info("%04d download skipped. already cached. (%s)", serialno, cached_path)
            return cached_path
        else:
            log.info("%04d no such content. downloading full file.", serialno)


        # download full file
        if not digest_type:
            digest_obj = hashlib.md5()
            digest_type = "md5"

        full_file = _do_download(tempdir, url, creds, logprefix="%04d.data " % (serialno,), digest_obj=digest_obj)
        computed_digest = digest_obj.hexdigest()

        if expected_digest:
            if computed_digest != expected_digest:
                log.error("%04d download digest %s:%s does not match expected %s:%s",
                          serialno, digest_type, computed_digest, digest_type, expected_digest)
                raise Exception("download failed. digest mismatch.")
            else:
                log.info("%04d digest %s:%s matches", serialno, digest_type, expected_digest)

        cached_path = cache.put(full_file, digest_type, computed_digest)
        return cached_path
    finally:
        shutil.rmtree(tempdir, ignore_errors=True)

def main():
    setup_logging()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outputdir", metavar="OUTPUTDIR", help="location where final files are placed")
    parser.add_argument("--threads", metavar="THREADS", type=int,
                        help="number of concurrent downloads")
    parser.add_argument("--workdir", metavar="WORKDIR", type=str, default="/tmp/",
                        help="working directory")
    parser.add_argument("--creds", metavar="CREDSFILE", type=str, default=None,
                        help="credentials file (yaml)")
    parser.add_argument("--keepgoing", help="continue after encountering errors", action="store_true",
                        default=False)

    args = parser.parse_args()

    if args.creds:
        creds_obj = load_creds(args.creds)
    else:
        creds_obj = {}

    creds_obj.setdefault('username', '')
    creds_obj.setdefault('password', '')

    try:
        if not os.path.exists(args.outputdir):
            os.makedirs(args.outputdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    cache = ContentCache(args.outputdir)

    if args.threads < 1:
        args.threads = 1

    total_done = 0
    total_submitted = 0
    errors = {}
    futures = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        for (lineno, line) in enumerate(sys.stdin):
            line = line.strip()
            if not line: continue
            toks = line.split()
            if len(toks) > 1:
                url, digest_type, digest_url = toks[0:3]
            else:
                url, digest_type, digest_url = toks[0], None, None

            fut = executor.submit(download_content, cache, args.workdir, url, digest_type, digest_url, creds_obj, lineno + 1)
            total_submitted += 1
            futures[fut] = (lineno + 1, url, digest_url)

        for fut in concurrent.futures.as_completed(futures):
            lineno, url, digest_url = futures[fut]
            total_done +=1
            log.debug("Download task #%04d terminated (%s)", lineno, url)
            log.info("Finished %d/%d tasks. (Total %5.2f%%)", total_done, total_submitted, (total_done *100.0)/total_submitted)
            try:
                cached_path = fut.result()
                log.info("Download task #%04d succeeded (%s)", lineno, url)
            except concurrent.futures.CancelledError:
                errors[lineno] = True
                log.info("Download task #%04d was cancelled (%s)", lineno, url)
                pass
            except Exception as exc:
                errors[lineno] = True
                log.error("Download %04d failed (%s)", lineno, url, exc_info=exc)
                if not args.keepgoing:
                    for x in futures:
                        x.cancel()

    if len(errors) > 0:
        raise Exception("%d downloads encountered errors." % (len(errors),))
    return 0

if __name__ == "__main__":
    main()
    sys.exit(0)
