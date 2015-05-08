#!/usr/bin/env python3

from Bio.SeqIO import parse, write
import sys

def main():
    for rec in parse(sys.argv[1], 'fastq'):
        write(rec, sys.stdout, 'fasta')

if __name__ == "__main__":
    main()
