#!/usr/bin/env python3


"""Process an input list (in a text file), containing paths or urls to fastq,
   bam, or bai files, and output a yaml-configuration file
   suitable for ingesting into the SNP calling pipeline Snakemake
   workflow. The program will extract sample name (and lane) information
   from basenames.

   Input files provided are expected to have one file path/or url per line,
   specifying either fastq (R1 and R2), md5, and bam/bai.

   The output is suitable as input to the snakemake snp-calling pipeline.

   e.g.:

     ./process_sample_list.py MY_LIST.TXT > samples.yaml


   use --help for complete usage.
"""

import sys
import re
import os.path
import argparse
import logging

from urllib.parse import urlparse

# HI.4530.001.index_7.ANN1317_R1.fastq.gz.md5  => HI.4530.001.index_7
LANE_RE = re.compile("^[hH][iI][.][0-9]+[.][0-9]+([.]([Ii]ndex|[Rr]ieseberg|[Bb]ioOHT|[Bb]ioO)_[0-9]+)?")

LOGGER = logging.getLogger()

def setupLogger():
    LOGGER.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

class Sample(object):
    """ each sample has the following attributes """
    def __init__(self, short_name, lane):
        super(Sample, self).__init__()
        self.short_name = short_name
        self.R1 = None
        self.R2 = None
        self.R1_MD5 = None
        self.R2_MD5 = None
        self.BAM = None
        self.BAI = None
        self.lane = lane

    def validate(self):
        if not self.short_name:
            raise ValueError("missing sample name")

        if (self.BAM and not self.BAI) or (self.BAI and not self.BAM):
            raise ValueError("sample %s has incomplete BAM/BAI pair" % (self.short_name))
        if self.BAM:
            return

        for k in ("R1" , "R2", "R1_MD5", "R2_MD5"):
            if not getattr(self, k):
                raise ValueError("sample %s is missing input file R1" % (self.short_name, k))

    def retrieve_bai(self, expanduser=True):
        if not self.BAM:
            return
        if self.BAI:
            return

        def _expand(pth):
            return os.path.expanduser(pth) if expanduser else pth

        tries = (
            self.BAM + ".bai",
            self.BAM[:-4] + ".bai"
        )
        for _try in tries:
            if os.path.exists(_expand(_try)):
                self.BAI = _try
                break

    def to_dict(self, ensure_scheme=False):
        def _ensure_scheme(pth):
            if not pth or not ensure_scheme:
                return pth
            if not pth.startswith("http:") and not pth.startswith("https:") and not pth.startswith("file:"):
                return "file://" + pth
            return pth

        return { "short_name": self.short_name or "",
                 "R1": _ensure_scheme(self.R1) or "",
                 "R2": _ensure_scheme(self.R2) or "",
                 "R1_MD5": _ensure_scheme(self.R1_MD5) or "",
                 "R2_MD5": _ensure_scheme(self.R2_MD5) or "",
                 "BAM": _ensure_scheme(self.BAM) or "",
                 "BAI": _ensure_scheme(self.BAI) or "",
                 "lane": self.lane
        }

def _identify_lane(fname):
    m = LANE_RE.match(fname)
    if m:
        return m.group(0)
    else:
        return None

def _identify_bam_bai(fname):
    """ extract the short name of the sample from its filename """
    basename = os.path.basename(fname)
    if not basename.endswith(".bai") and not basename.endswith(".bam"):
        return (None, None, None)

    lane = _identify_lane(basename)

    if basename.endswith(".bam.bai"):
        typ = "BAI"
        front, _ = basename.rsplit(".bam.bai", 1)
    elif basename.endswith(".bai"):
        typ = "BAI"
        front, _ = basename.rsplit(".bai", 1)
    else:
        typ = "BAM"
        front, _ = basename.rsplit(".bam", 1)
    leading_bits, short_name = front.rsplit(".", 1)
    return (lane, short_name, typ)

def identify_sample(fname):
    """ extracts the short name of the sample from the filename
        filename can be a URL or just a path on disk

    >>> identify_sample("https://genomequebec.mcgill.ca/nanuqMPS/readSetMd5Download/id/417446/type/READ_SET_FASTQ/filename/HI.4530.001.index_7.ANN1317_R1.fastq.gz.md5")
    ('HI.4530.001', 'ANN1317', 'R1_MD5')

    >>> identify_sample("https://genomequebec.mcgill.ca/nanuqMPS/readSetMd5Download/id/417445/type/READ_SET_FASTQ_PE/filename/HI.4530.001.Index_3.HK24_16_01_R_010_R2.fastq.gz.md5")
    ('HI.4530.001', 'HK24_16_01_R_010', 'R2_MD5')

    >>> identify_sample("https://genomequebec.mcgill.ca/nanuqMPS/fileDownload/id/417446/type/READ_SET_FASTQ/filename/HI.4530.001.index_7.ANN1317_R1.fastq.gz")
    ('HI.4530.001', 'ANN1317', 'R1')

    >>> identify_sample("https://genomequebec.mcgill.ca/nanuqMPS/fileDownload/id/417449/type/READ_SET_FASTQ_PE/filename/HI.4530.001.Index_16.PET0520_R2.fastq.gz")
    ('HI.4530.001', 'PET0520', 'R2')
    """

    parsed = urlparse(fname, allow_fragments=True)
    pathname = parsed.path
    basename = os.path.basename(pathname)

    supported_ext = (
        ".fastq", ".fastq.md5",
        ".fastq.gz", ".fastq.gz.md5",
        ".bam", ".bai"
    )
    for suffix in supported_ext:
        if basename.endswith(suffix):
            break
    else:
        return (None, None, None)

    if fname.endswith(".bam") or fname.endswith(".bai"):
        return _identify_bam_bai(basename)

    lane = _identify_lane(basename)

    front, tail = basename.rsplit(".fastq.", 1)
    leading_bits, short_name_RX = front.rsplit(".", 1)


    if not (short_name_RX.endswith("_R2") or short_name_RX.endswith("_R1")):
        raise Exception("cannot extract short sample name from url `%s`. Expecting to find _R1 or _R2." %
                        (fname,))

    short_name, rtype = short_name_RX.rsplit("_", 1)

    if tail.endswith(".md5"):
        # e.g. ("PET0520", "R1_MD5")
        return (lane, short_name, rtype.upper() + "_" + "MD5")
    else:
        # e.g. ("PET0520", "R1")
        return (lane, short_name, rtype.upper())


def _output_yaml(all_samples, sample_prefix, problems=None, outfd=sys.stdout):
    """output samples as YAML"""

    if problems:
        outfd.write("#Warning. %d lines were skipped.\n" % (len(problems),))
        for prob in problems:
            outfd.write("#L%04d %s\n" % (prob[0], prob[1]))
        outfd.write("#\n")

    outfd.write("samples:\n")
    for sample_key in sorted(all_samples.keys()):
        sample_name, lane = sample_key
        sample = all_samples[sample_key]
        vals = sample.to_dict(ensure_scheme=True)
        vals["short_name"] = sample_prefix + vals["short_name"]
        outfd.write(    "  - name: %(short_name)s\n" % vals)
        if 'lane' in vals:
            outfd.write("    lane: %(lane)s\n" % vals)
        outfd.write(    "    locations:\n")
        if sample.R1 or sample.R2:
            outfd.write("      fastq:\n")
            outfd.write("        - %(R1)s\n" % vals)
            outfd.write("        - %(R2)s\n" % vals)
            outfd.write("      md5:\n")
            outfd.write("        - %(R1_MD5)s\n" % vals)
            outfd.write("        - %(R2_MD5)s\n" % vals)
        if sample.BAM:
            outfd.write("      bam: %(BAM)s\n" % vals)
            outfd.write("      bai: %(BAI)s\n" % vals)

def _output_txt(all_samples, sample_prefix, problems=None, outfd=sys.stdout):
    """output samples in a form that makes it easy to download files:

       URL md5 MD5_URL

    """
    if problems:
        outfd.write("#Warning. %d lines were skipped.\n" % (len(problems),))
        for prob in problems:
            outfd.write("#L%04d %s\n" % (prob[0], prob[1]))
        outfd.write("#\n")

    for sample_name in sorted(all_samples.keys()):
        outfd.write("%(R1)s md5 %(R1_MD5)s\n" % all_samples[sample_name].to_dict())
        outfd.write("%(R2)s md5 %(R2_MD5)s\n" % all_samples[sample_name].to_dict())

def process_input_file(fpath, sample_prefix, output_fn, addbai=False, expanduser=True, **args):
    all_samples = {}

    problem_lines = []
    with open(fpath, "r") as fd:
        for lineno, line in enumerate(fd):
            if lineno % 100 == 0 and lineno > 0:
                LOGGER.debug("processed %d lines...", lineno)

            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                lane, sample_name, file_type = identify_sample(line)
                if sample_name is None:
                    continue
            except Exception as exc:
                LOGGER.warn("Problem. Skipping line %d: %s" % (lineno + 1, line))
                problem_lines.append((lineno, line))
                continue
            # we assume each sample name will have at most one bam/bai
            # and at most one fastq read set per lane.
            sample_key = (sample_name, lane)
            sample = all_samples.setdefault(sample_key, Sample(sample_name, lane))

            prev = getattr(sample, file_type)
            if prev is not None:
                raise Exception("Sample %s lane %s, has duplicate values for entry %s" % (
                    sample_name, lane, file_type))
            setattr(sample, file_type, line)

    LOGGER.info("validating %d samples...", len(all_samples))

    for sampleno, (sample_key, sample) in enumerate(all_samples.items()):
        if sampleno % 100 == 0 and sampleno > 0:
            LOGGER.debug("validated %d samples...", sampleno)

        if addbai and sample.BAM and not sample.BAI:
            sample.retrieve_bai(expanduser=expanduser)
        sample.validate()
    LOGGER.info("all %d samples validated.", len(all_samples))

    return output_fn(all_samples, sample_prefix, problem_lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("INPUTFILE", help="the list of input samples (fastq(.gz)) in URLs or paths")
    parser.add_argument("--mode", help="output mode", type=str, choices=["yaml", "urls"], default="yaml")
    parser.add_argument("--nameprefix", help="add prefix to sample short name", type=str, default="")
    parser.add_argument("--addbai", help="if the list provides only a bam, search for a bai file based on the bam name.",
                        action="store_true", default=False)
    parser.add_argument("--no-expanduser", help="when checking for existence of files, do not expand ~ and ~user inside paths",
                        dest="expanduser", action="store_false", default=True)
    args = parser.parse_args()
    setupLogger()

    mode = args.mode
    if mode == "yaml":
        process_input_file(args.INPUTFILE, args.nameprefix, _output_yaml, **vars(args))
    elif mode == "urls":
        process_input_file(args.INPUTFILE, args.nameprefix, _output_txt, **vars(args))
    else:
        raise Exception("invalid mode")

    sys.exit(0)
