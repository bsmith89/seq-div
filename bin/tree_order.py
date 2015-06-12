#!/usr/bin/env python3

from Bio.Phylo import read
from Bio.SeqIO import index, write
import sys


def main():
    tree = read(sys.argv[1], 'newick')
    seqs = index(sys.argv[2], 'fasta')
    if not tree.rooted:
        tree.root_at_midpoint()
    tree.ladderize(reverse=True)
    for leaf in tree.get_terminals():
        write(seqs[leaf.name], sys.stdout, 'fasta')

if __name__ == "__main__":
    main()
