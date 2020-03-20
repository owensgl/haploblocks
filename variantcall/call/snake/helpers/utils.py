
import os, os.path
import errno
import tempfile
import hashlib
import io
import traceback

class Cas(object):
    """ Content Addressible Storage.

        Storage device/directory where all files are stored under their content address (hash).

        put(file) -> key
        get(key)  -> file

        two files with the same content will have the same key
    """
    def __init__(self, cachedir):
        self.cachedir = cachedir

    def _key(self, digest_type, hexdigest):
        if digest_type == "sha1":
            return hexdigest.lower()
        return digest_type.upper() + "_" + hexdigest.lower()

    def _cache_entry(self, digest_type, digest):
        return os.path.join(self.cachedir, digest_type.strip().upper() + "_" + digest.strip().lower())

    def _get(self, digest_type, hexdigest):
        return open(self._cache_entry(digest_type, hexdigest), "r")

    def get(self, name):
        """gets an open file based on key. raises erorr if file is not found."""
        toks = name.split("_", 1)
        if len(toks) == 1:
            digest_type, hexdigest = "sha1", toks[0]
        else:
            digest_type, hexdigest = toks
        return self._get(digest_type, hexdigest)

    def put_text(self, text, digest_type='sha1'):
        """stores text string. returns a key that can be used to retrieve content later"""

        with io.BytesIO(text.encode('utf-8')) as datafile:
            return self._put(datafile, digest_type=digest_type)

    def put_stream(self, stream, digest_type='sha1'):
        """stores file-like object stream.
           returns a key that can be used to retrieve content later"""
        return self._put(stream, digest_type=digest_type)

    def _put(self, stream, digest_type):
        digest_obj = getattr(hashlib, digest_type)()

        with tempfile.NamedTemporaryFile(prefix="tmp.digest.", dir=self.cachedir, delete=False) as digest_fd:
            chunk = stream.read(8*1024)
            while chunk:
                digest_obj.update(chunk)
                digest_fd.write(chunk)
                chunk = stream.read(8*1024)

            hexdigest = digest_obj.hexdigest()
            dst_name = self._cache_entry(digest_type, hexdigest)

        # avoids modifying timestamps
        if not os.path.exists(dst_name):
            os.rename(digest_fd.name, dst_name)
        else:
            os.remove(digest_fd.name)

        return self._key(digest_type, hexdigest)

def log_exc(f, logger):
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as exc:
            trace = traceback.format_exc(10)
            logger.error("%s: %s" % (exc, trace))
            raise
    return wrapper

def mkdirp(pathname, mode=0o755):
    """mkdir -p <pathname>"""
    try: os.makedirs(pathname, mode)
    except OSError as ose:
        if ose.errno != errno.EEXIST:
            pass

