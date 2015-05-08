#!/usr/bin/env python3
"""Construct FASTQ files from PHRED output."""

from Bio.SeqIO import read, write
import sys

def main():
    seq = read(sys.argv[1], 'fasta')
    qual = read(sys.argv[2], 'qual')
    seq.letter_annotations = qual.letter_annotations
    write(seq, sys.stdout, 'fastq')


if __name__ == "__main__":
    main()
