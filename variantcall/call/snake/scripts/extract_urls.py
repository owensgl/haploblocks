#!/bin/env python3

#
# extracts urls from a sample file
#



import yaml
import sys

def do_file(fd):
    obj = yaml.load(fd)

    for sample in obj['samples']:
        if not "locations" in sample: continue
        fastqs = sample["locations"].get("fastq", [])
        md5s = sample["locations"].get("md5", [])
        for i, fastq in enumerate(fastqs):
            md5f = md5s[i] if len(md5s) > i else None
            if md5f:
                yield (fastq, md5f)
            else:
                yield (fastq, None)


def output_record(fastq, md5):
    if md5:
        print("%s md5 %s" % (fastq, md5))
    else:
        print("%s" % (fastq,))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        for fastq, md5 in do_file(sys.stdin):
            output_record(fastq, md5)
    else:
        for fname in sys.argv[1:]:
            with open(fname, "r") as fd:
                for fastq, md5 in do_file(fd):
                    output_record(fastq, md5)


