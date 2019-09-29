#!/usr/bin/env python3
#
# Counts the number of basepairs and intervals defined in the BED file
# passed via stdin
#


import sys
from collections import OrderedDict

if __name__ == "__main__":
    db = OrderedDict()
    intervals = 0
    totalbp = 0
    for line in sys.stdin:
        toks = line.strip().split()
        chrom, start, end = toks
        start = int(start, 10)
        end = int(end, 10)
        db[chrom] = db.get(chrom, 0) + (end - start)
        totalbp += end - start
        intervals += 1
    for chrom in db:
        print("%s %d" % (chrom, db[chrom]))
    print("num_intervals %d" % (intervals,))
    print("total_bp %d" % (totalbp,))

        
        
        
