#!/usr/bin/env python3
"""Parse primersearch output files into tables.

TODO: Extract information about the sequence where the primers hit."""

import sys
from Bio.Emboss.PrimerSearch import read
from pandas import DataFrame
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna


def psearch_primer_to_iupac(string):
    pieces = []
    for piece in string.split("]"):
        if "[" in piece:
            unambig, ambig = piece.split("[")
            pieces.append(unambig)
            pieces.append(ambig[-1])
        else:
            pieces.append(piece)
    out1 = ''.join(pieces)
    out2 = 'N'.join(out1.split('?'))
    return Seq(out2, alphabet=ambiguous_dna)

def parse_primer_sets(handle):
    pairs = {}
    for line in handle:
        name, forward, reverse = tuple(word.strip() for word in line.split())
        pairs[name] = (SeqRecord(id="forward",
                           seq=Seq(forward.upper(), alphabet=ambiguous_dna)),
                       SeqRecord(id="reverse",
                           seq=Seq(reverse.upper(), alphabet=ambiguous_dna)))
    return pairs

def extract_info(amplicon, primers):
    name = amplicon.seq_id
    length = amplicon.length
    primer_start, primer_stop = None, None
    for primer in amplicon.binding_sites:
        index, mismatches = amplicon.binding_sites[primer]
        if index[0] != '[':  # PrimerSearch notation for the complementary strand
            start, mis_start, primer_start = int(index) - 1, mismatches, primer
        elif index[0] == '[':
            stop, mis_stop, primer_stop = int(index[1:-1]), mismatches, primer
    if not primer_start and primer_stop:
        print(primer)
    primer_start = psearch_primer_to_iupac(primer_start)
    primer_stop = psearch_primer_to_iupac(primer_stop)
    primer_start_len = len(primer_start)
    primer_stop_len = len(primer_stop)
    if (primer_start == primers[0].seq) and (primer_stop == primers[1].seq):
        primer_start_name = primers[0].id
        primer_stop_name = primers[1].id
    elif (primer_stop == primers[0].seq) and (primer_start == primers[1].seq):
        primer_stop_name = primers[0].id
        primer_start_name = primers[1].id
    else:
        raise ValueError("Primers {} do not match primers {}".\
                             format((primer_start, primer_stop), (primers)))
    return (name,
            start, stop,
            mis_start, mis_stop,
            primer_start_name, primer_stop_name,
            primer_start_len, primer_stop_len)

def mung_file(psearch_handle, primers_handle):
    records = []
    primer_sets = parse_primer_sets(primers_handle)
    amplifiers = read(psearch_handle).amplifiers
    for primer_set in amplifiers:
        for amplicon in amplifiers[primer_set]:
            records.append([primer_set] + \
                list(extract_info(amplicon, primer_sets[primer_set])))
    return DataFrame(records, columns=('primer_set', 'template',
                                       'start', 'stop',
                                       'mis_start', 'mis_stop',
                                       'primer_start', 'primer_stop',
                                       'len_start', 'len_stop'))

def main():
    with open(sys.argv[1]) as primers_handle, \
         open(sys.argv[2]) as psearch_handle:
        data = mung_file(psearch_handle, primers_handle)
        data.to_csv(sys.stdout, sep='\t', index=False)

if __name__ == '__main__':
    main()
