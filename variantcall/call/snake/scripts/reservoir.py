#!/usr/bin/env python3

"""
   Sample N random lines from the input in such a way that
   each line from the input is selected with equal probability.

   Input feed is stdin, output is stdout.

   Lines are output only once EOF is reached.
"""

import os
import sys
import argparse
import random
import time

class Unbuffered(object):
    """saves us from calling python -u"""
    def __init__(self, stream):
       self.stream = stream
    def write(self, data):
       if self.stream.write(data):
           self.stream.flush()
    def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
    def __getattr__(self, attr):
       return getattr(self.stream, attr)

class DataException(Exception):
    pass

def reservoir(infd, outfd, progressfd, keep_count, keep_cols=None, sep=None, has_header=True, output_sep="\t"):
    if keep_count < 1:
        return

    selected = []
    col_ids = None
    line = None

    class progress_ctx: pass

    progress_ctx.start_ts = time.time()
    progress_ctx.last_ts = progress_ctx.start_ts
    progress_ctx.last_lineno = 0

    def _show_progress(ctx, lineno, line):
        now = time.time()
        if (now - ctx.last_ts < 5.0):
            return
        avg_speed = (lineno - ctx.last_lineno) / (now - ctx.last_ts + 0.000001)
        progressfd.write("[%8.3f] reservoir progress: line %10d: %20s... (%9.3f lines/min)\n" % (
            (now - ctx.start_ts), lineno , _line_data(line).strip()[0:20], avg_speed * 60.0))
        ctx.last_ts = now
        ctx.last_lineno = lineno

    def _line_data(line):
        if col_ids:
            toks = line.split(sep=sep)
            try:
                return output_sep.join(tuple(toks[idx] for idx in col_ids)) + "\n"
            except IndexError as ie:
                raise DataException("Missing columns in line {!r}".format(line))
        else:
            return line

    if not has_header:
        if keep_cols:
            # numeric column names only
            col_ids = [int(colname, 10) - 1 for colname in keep_cols]
        out_header = None
    else:
        for line in infd: break # one line
        if not line:
            raise DataException("Expected a header. Got EOF.")
        if keep_cols:
            colmap = {colname: colidx for (colidx, colname) in enumerate(line.split(sep=sep))}
            try:
                col_ids = [(colmap[colname] if colname[0] not in "1234567890" else int(colname, 10) - 1) for colname in keep_cols]
            except KeyError as ke:
                raise DataException("column name not found: {}".format(ke))
        out_header = _line_data(line)

    for lineno, line in enumerate(infd):
        if lineno % 100000 == 0:
            _show_progress(progress_ctx, lineno, line)

        if len(selected) < keep_count:
            selected.append(_line_data(line))
            continue

        # with probability keep_count/i, keep item
        random_index = random.randint(0, lineno)
        if random_index < keep_count:
            selected[random_index] = _line_data(line)

    if out_header:
        outfd.write(out_header)
    for datum in selected:
        outfd.write(datum)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("count", metavar="COUNT", type=int,
                        help="the number of records to keep")
    parser.add_argument("--col", metavar="COL", type=str, action='append',
                        help="add column with header COL to the set of kept columns. COL is the column name or a 1-based numeric value.")
    parser.add_argument("--sep", metavar="SEP", type=str, default=None,
                        help="when selecting fields, use SEP as the separator. default is any whitespace character")
    parser.add_argument("--no-header", dest="header", action="store_false", default=True,
                        help="set this if data starts on the first line. the default is to assume a header")
    parser.add_argument("--osep", metavar="OSEP", type=str, default="\t",
                        help="output field separator")
    args = parser.parse_args()

    if args.count < 0:
        count = 0

    random.seed()
    try:
        unbuffered_err = Unbuffered(sys.stderr)
        reservoir(sys.stdin, sys.stdout, unbuffered_err, args.count, args.col, sep=args.sep, has_header=args.header, output_sep=args.osep)
    except DataException as de:
        sys.stderr.write("reservoir.py: error: {}\n".format(de))
        sys.exit(1)

    sys.exit(0)

