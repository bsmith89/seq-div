#!/usr/bin/env python3
"""Take only good sequences and trim them down to an internal primer set.

Takes 3 command line arguments:
    (1) primersearch hits table as output from bin/mung_primersearch.py
    (2) sequence file in fasta
    (3) the first n in-frame letters of the correctly oriented sequence
        (probably the first n letters of the search primer that you'd like to
         be at the 5' end.)

Defaults to 3 mismatches per primer allowed.

TODO: Fix up this script a LOT.


"""

import sys
import pandas as pd
from Bio.SeqIO import index, write
from warnings import warn

def main():
    hits_table = pd.read_table(sys.argv[1])
    hits_table.template = hits_table.template.astype(str)
    good_hits = hits_table[(hits_table.mismatch_start <= 2) &
                           (hits_table.mismatch_stop <= 2)]
    seq_index = index(sys.argv[2], 'fasta')
    out = []
    for _, hit in good_hits.iterrows():
        rec = seq_index[hit['template']]
        seq = rec.seq[(hit['start']):(hit['stop'])]
        start_seq = sys.argv[3]
        # Re-orient sequences based on argv[3] sequence.
        if seq[:len(start_seq)].upper() == start_seq.upper():
            rec.seq = seq
        else:
            rec.seq = seq.reverse_complement()
            if rec.seq[:len(start_seq)].upper() != start_seq.upper():
                warn("Correct orientation of {} could not be found. "
                     "Dropping sequence.".format(rec.id))
                continue
        # TODO: Remove primers
        # Set frame
        while len(rec.seq) % 3 != 0:
            rec.seq = rec.seq[:-1]
        out.append(rec)
    write(out, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
