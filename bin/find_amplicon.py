#!/usr/bin/env python3
"""Trim reads based on where a primer hits."""


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

logger = logging.getLogger(__name__)
logging.captureWarnings(True)

def parse_primer_sets(handle):
    pairs = {}
    for line in handle:
        name, forward, reverse = tuple(word.strip() for word in line.split())
        pairs[name] = (SeqRecord(id=name + "_forward",
                           seq=Seq(forward, alphabet=IUPACAmbiguousDNA)),
                       SeqRecord(id=name + "_reverse",
                           seq=Seq(reverse, alphabet=IUPACAmbiguousDNA)))
    return pairs

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
    if len(hit_info) == 0:
        raise NoHitsError(template_id)
    elif len(hit_info) > 1:
        if not best:
            raise MultipleHitsError(template_id)
        else:
            hit_info['mismatch_sum'] = hit_info.mismatch_start + hit_info.mismatch_stop
            hit_info = hit_info.sort('mismatch_sum').iloc[0:1]
    else:
        pass  # hit_info is already the one hit
    return hit_info

def _get_extra_args():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*cli.POS_GROUP)
    g.add_argument("primers_handle", metavar='PRIMERS',
                   type=argparse.FileType('r'),
                   help="primers file")
    g.add_argument("table_handle", metavar='AMPLICONS',
                   type=argparse.FileType('r'),
                   help="parsed primersearch results")
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

def main():
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    primer_sets = parse_primer_sets(args.primers_handle)
    hits_table = read_table(args.table_handle)
    recs = parse(args.in_handle, args.fmt_infile)
    for rec in recs:
        try:
            hit_info = get_hit_info(rec.id, hits_table)
        except (NoHitsError, MultipleHitsError) as err:
            warn(str(err))
            rec = rec[0:0]
            write(rec, args.out_handle, args.fmt_outfile)
            continue
        hit_info = {col:hit_info[col].iloc[0] for col in hit_info.columns}
        hit_primer_set = hit_info['primer_set']
        fwd_primer = str(primer_sets[hit_primer_set][0].seq).upper()
        rev_primer = str(primer_sets[hit_primer_set][1].seq).upper()
        hit_primer_start = hit_info['primer_start'].upper()
        hit_primer_stop = hit_info['primer_stop'].upper()
        hit_start = hit_info['start']
        hit_stop = hit_info['stop']
        trim_start = hit_start + len(fwd_primer)
        trim_stop = hit_stop - len(rev_primer)
        if hit_primer_start != fwd_primer:
            assert hit_primer_stop == fwd_primer
            out_rec = rec[trim_start:trim_stop].reverse_complement()
            out_rec.id = rec.id
            out_rec.description = rec.description
        else:
            pass
            out_rec = rec[trim_start:trim_stop]
        write(out_rec, args.out_handle, args.fmt_outfile)

if __name__ == "__main__":
    main()
