#!/usr/bin/env python3
"""Remove sequences of length other than the modal length."""

from collections import Counter, OrderedDict
from Bio.SeqIO import parse, write
from utils.lib import cli
import sys
import argparse
import logging
from warnings import warn
from copy import deepcopy

logger = logging.getLogger(__name__)
logging.captureWarnings(True)

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              cli.get_seq_in_parser(),
                                              ])
    args = parser.parse_args(argv[1:])
    return args

def main():
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    all_recs = OrderedDict()
    for rec in parse(args.in_handle, args.fmt_infile):
        all_recs[rec] = len(rec)

    mode = Counter(all_recs.values()).most_common()[0][0]
    for rec in all_recs:
        if all_recs[rec] == mode:
            write(rec, args.out_handle, args.fmt_outfile)
        else:
            warn(cli.DropSequenceWarning(
                "{} had length {}, not {}".format(rec.id, len(rec), mode)))

if __name__ == '__main__':
    main()
