#!/usr/bin/env python3
"""Trim sequences to a high quality center.

"""

import sys
from Bio.SeqIO import parse, write
from statistics import mean
from utils.lib import cli
import argparse
import logging

DEFAULT_QUAL_THRESH = 40

logger = logging.getLogger(__name__)

def quality_trim(rec, end, avg):
    start = 0
    stop = len(rec)
    quality = rec.letter_annotations['phred_quality']
    while (quality[stop - 1] < end) and (stop > 0):
        stop -= 1
    while (quality[start] < end) and (stop - start > 0):
        start += 1
    while (stop - start > 0) and (mean(quality[start:stop]) < avg):
        if quality[stop - 1] < quality[start]:
            stop -= 1
        else:  # In the case of a tie, the 5' end is truncated.
            start += 1
    return rec[start:stop]

def _get_quality_args_parser():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*cli.PAR_GROUP)
    g.add_argument('--quality-threshold', '-q', type=int,
                   default=DEFAULT_QUAL_THRESH,
                   help=("mean PHRED-quality threshold "
                         "DEFAULT: {}").format(DEFAULT_QUAL_THRESH))
    return p

def _get_fastq_in_parser():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*cli.POS_GROUP)
    g.add_argument("in_handle",
                   type=argparse.FileType('r'),
                   metavar="SEQUENCE", default=sys.stdin,
                   help=("sequence file"))
    return p


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              _get_quality_args_parser(),
                                              _get_fastq_in_parser(),
                                              ])
    args = parser.parse_args(argv[1:])
    return args

def main():
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for rec in parse(args.in_handle, 'fastq'):
        write(quality_trim(rec,
                           args.quality_threshold,
                           args.quality_threshold),
              args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
