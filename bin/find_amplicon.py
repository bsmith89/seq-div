#!/usr/bin/env python3
"""Trim reads based on where a primer hits.

TODO: Calculate the temperature of melting for each hit and determine
if it would be amplified.

"""


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqIO import parse, write
import sys
from pandas import read_table, Series
from warnings import warn
from utils.lib import cli
import argparse
import logging
from copy import copy


logger = logging.getLogger(__name__)
logging.captureWarnings(True)

class Error(Exception):
    pass


class NoHitsError(Error):
    def __init__(self, template_id):
        self.template_id = template_id
    def __str__(self):
        return "Hits to {} not found.".format(self.template_id)


class MultipleHitsError(Error):
    def __init__(self, template_id):
        self.template_id = template_id
    def __str__(self):
        return "Multiple hits to {} found.".format(self.template_id)


def get_hit_info(template_id, hits_table, best=True):
    hit_info = hits_table[hits_table.template == template_id]
    if (len(hit_info) > 1) and (not best):
        raise NoHitsError(template_id)
    elif (len(hit_info) == 0):
        raise NoHitsError(template_id)
    else:
        return hit_info.ix[hit_info.sort('mis_sum').index[0]]

def _get_extra_args():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*cli.POS_GROUP)
    g.add_argument("table_handle", metavar='PSEARCH',
                   type=argparse.FileType('r'),
                   help="parsed primersearch results")
    h = p.add_argument_group(*cli.PAR_GROUP)
    h.add_argument("--max-mismatch", type=int,
                   default=None,
                   help=("maximum mismatches allowed. "
                         "DEFAULT: no limit"))
    h.add_argument("--primer-set", type=str,
                   default=None,
                   help=("name of the primer set to use. DEFAULT: no consideration"))
    h.add_argument("--trim-primers", action="store_true",
                    help="include sequence region covered by primer hit.")
    return p


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              _get_extra_args(),
                                              cli.get_seq_in_parser()
                                              ])
    args = parser.parse_args(argv[1:])
    return args

def get_amplicon(rec, hits, trim_primers=False):
    out_rec = copy(rec)
    try:
        hit_info = get_hit_info(rec.id, hits)
    except NoHitsError as err:
        out_rec = out_rec[:0]
        return out_rec, None
    if trim_primers:
        trim_start = hit_info['start'] + hit_info['len_start']
        trim_stop = hit_info['stop'] + hit_info['len_stop'] - 1
    else:
        trim_start = hit_info['start']
        trim_stop = hit_info['stop'] - 1
    out_rec = out_rec[trim_start:len(out_rec) - trim_stop]
    if hit_info['primer_start'] == 'reverse':
        id, description = out_rec.id, out_rec.description
        out_rec = out_rec.reverse_complement()
        out_rec.id = id
        out_rec.description = description
    return out_rec, hit_info

def main():
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    hits = read_table(args.table_handle)
    hits['mis_sum'] = hits.mis_start + hits.mis_stop
    if args.max_mismatch:
        hits = hits[hits.mis_sum <= args.max_mismatch]
    if args.primer_set:
        hits = hits[hits.primer_set == args.primer_set]

    recs = parse(args.in_handle, args.fmt_infile)
    for rec in recs:
        amplicon, hit_info = get_amplicon(rec, hits,
                                          trim_primers=args.trim_primers)
        logger.debug(hit_info)
        if (type(hit_info) == type(None)) and args.drop:
            warn(cli.DropSequenceWarning("No hit found for {rec.id}".format(rec=rec)))
        else:
            write(amplicon, args.out_handle, args.fmt_outfile)

if __name__ == "__main__":
    main()
