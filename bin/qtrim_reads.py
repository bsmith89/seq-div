#!/usr/bin/env python3
"""Trim sequences to a high center.  Throw them away entirely if they're
too short.


"""

import sys
from Bio.SeqIO import parse, write
from statistics import mean

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

def main():
    for rec in parse(sys.argv[1], 'fastq'):
        write(quality_trim(rec, 40, 40), sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
