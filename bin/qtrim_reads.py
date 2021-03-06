#!/usr/bin/env python3
"""Trim sequences to a high quality center.

"""

import sys
from Bio.SeqIO import parse, write
from statistics import mean
from utils.lib import cli
import argparse
import logging
from warnings import warn
from copy import deepcopy

DEFAULT_QUAL_THRESH = 40
DEFAULT_MIN_LEN = 9

logger = logging.getLogger(__name__)
logging.captureWarnings(True)

def quality_trim(rec, thr, just_filter=False, keep_columns=False):
    if keep_columns:
        raise NotImplementedError("This feature is not yet implemented.")

    rec = deepcopy(rec)
    if len(rec.seq) == 0:
        return rec

    start = 0
    stop = len(rec)
    quality = rec.letter_annotations['phred_quality']
    if just_filter and min(quality) < thr:
        stop = 0  # Trim the whole thing.


    while (quality[stop - 1] < thr) and (stop > 0):
        stop -= 1
    while (quality[start] < thr) and (stop - start > 0):
        start += 1
    while (stop - start > 0) and (mean(quality[start:stop]) < thr):
        if quality[stop - 1] < quality[start]:
            stop -= 1
        else:  # In the case of a tie, the 5' thr is truncated.
            start += 1
    return rec[start:stop]

def _get_quality_args_parser():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*cli.PAR_GROUP)
    g.add_argument('--quality-threshold', '-q', type=int,
                   default=DEFAULT_QUAL_THRESH,
                   help=("PHRED-quality threshold "
                         "DEFAULT: {}").format(DEFAULT_QUAL_THRESH))
    g.add_argument('--min-length', '-l', type=int,
                   default=DEFAULT_MIN_LEN,
                   help=("minimum length sequence; "
                         "makes noise only unless --drop "
                         "DEFAULT: {}").format(DEFAULT_MIN_LEN))
    g.add_argument('--keep-columns', action='store_true',
                   help=("replace trimmed nucleotides with '-'; "
                         "this makes sense if the reads are already aligned, "
                         "because they have been trimmed based on primers, "
                         "and no InDels are expected within the amplicon."))
    g.add_argument('--just-filter', action='store_true',
                   help=("drop sequences with mean quality less than "
                         "QUALITY_THRESHOLD, do not attempt to trim."))
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

    for rec_in in parse(args.in_handle, 'fastq'):
        logger.debug(rec_in)
        rec_out = quality_trim(rec_in, args.quality_threshold,
                               keep_columns=args.keep_columns)
        length = len(rec_out.seq)
        if length < args.min_length:
            warn(("Length of sequence {} less than threshold. "
                  "{} < {}. Dropping.").\
                     format(rec_out.id, length, args.min_length),
                 cli.DropSequenceWarning)
        else:
            write(rec_out, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
