#!/usr/bin/env python3

from Bio.SeqIO import parse, write
import sys
from collections import Counter

def in_frame(seq, frame):
    if frame < 0:
        seq = seq.reverse_complement()
        seq = seq[((-frame) - 1):]
    else:
        seq = seq[(frame - 1):]
    while len(seq) % 3 != 0:
        seq = seq[:-1]
    return seq

def inframe(seq):
    candidates = []
    for frame in [-3, -2, -1, 1, 2, 3]:
        candidate = in_frame(seq, frame)
        candidates.append((Counter(candidate.translate())['*'], candidate))
    stops, result = sorted(candidates)[0]
    # sys.stderr.write("{}\n".format(result.translate()))
    return result



def main():
    for rec in parse(sys.argv[1], 'fasta'):
        rec.seq = inframe(rec.seq)
        write(rec, sys.stdout, 'fasta')

if __name__ == "__main__":
    main()
