#!/usr/bin/env python3

"""
   window-creator-all.py [--size 1000000] [REFERENCEFILE | - ]

   reads a reference file (fa), to obtain chromosome headers and their
   sizes, and prints NAME:START-END where start and end are basepair
   offsets _within_ chromosome NAME.

   You may accelerate the output by using grep as a pre-filtering step:

      grep ">" referencefile | window-creator-all.py -

   this script is pretty much the same as window_creator.pl, except
   that it is meant to run a single time, and spits out information
   in space separated values.

   outputs contigs to stdout. one window per line:

   CHROMi STARTi ENDi CHROMLENi
   CHROMj STARTj ENDj CHROMLENj
   ...

   the last column, is a convenience to retrieve the length of the entire
   chromosome.
"""

import sys
import argparse
import yaml

def window_split(chrom, chromlen, size):
    num_windows = chromlen // size
    rem = chromlen % size
    for i in range(0, num_windows):
        yield {
            'chrom': chrom,
            'start': i*size + 1,
            'end': (i+1)*size,
            'clen': chromlen
        }
    if rem > 0:
        yield {
            'chrom': chrom,
            'start': 1+(num_windows*size),
            'end': chromlen,
            'clen': chromlen
        }

def chromosome_headers(infd, is_index):
    """ generator of (chromname, length_bp) tuples over input file """

    def bork(lineno, msg):
        raise Exception("Error occurred at Line %d: %s" % (lineno, msg))

    if is_index:
        # every line in a fasta index is a valid header
        for i, line in enumerate(infd):
            if not line: continue
            tokens = line.strip().split()
            if len(tokens) < 3:
                bork(i+1, "expected fasta index to have at least three columns")

            chrom = tokens[0].strip()

            try:
                chromlen = int(tokens[1], 10)
            except ValueError:
                bork(i+1, "invalid length column")
            yield (chrom, chromlen)
    else:
        for i, line in enumerate(infd):
            if not line.startswith(">"):
                continue
            tokens = line.strip().split()
            if len(tokens) < 2:
                bork(i+1, "expected fasta file to have annotations with each sequence name")
            chrom = tokens[0][1:]
            chromlen = tokens[1]
            if not chromlen.startswith("len="):
                bork(i+1, "expected fasta file to have len=N annotations with each sequence name")
            chromlen = chromlen[len("len="):]
            try:
                chromlen = int(chromlen, 10)
            except ValueError:
                bork(i+1, "invalid number found in len=N annotation")
            yield (chrom, chromlen)

def main(infd, outfd, size=1000000, is_index=False, **_):
    for chrom, chromlen in chromosome_headers(infd, is_index=is_index):
        for window in window_split(chrom, chromlen, size):
            outfd.write("%(chrom)s %(start)d %(end)d %(clen)d\n" % window)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("reference", metavar="REFERENCE", default="-",
                        help="the path to the reference (.fa). defaults to '-'"
                        " (for reading from stdin)")
    parser.add_argument("--size", metavar="SIZE", type=int,
                        help="size of windows created. in basepairs",
                        default=1000000)
    parser.add_argument("--faidx", action="store_true", dest='is_index',
                        help="set to true if the input is the fasta index")

    args = parser.parse_args()

    ref = args.reference
    if ref == "-":
        fd = sys.stdin
    else:
        fd = open(args.reference, "r")
    try:
        main(fd, sys.stdout, **vars(args))
    finally:
        fd.close()

    sys.exit(0)

