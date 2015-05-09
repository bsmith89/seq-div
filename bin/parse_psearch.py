#!/usr/bin/env python3
"Parse primersearch output files into tables."

import sys
from Bio.Emboss.PrimerSearch import read
from pandas import DataFrame

def psearch_primer_to_iupac(string):
    out = []
    for piece in string.split("]"):
        if "[" in piece:
            unambig, ambig = piece.split("[")
            out.append(unambig)
            out.append(ambig[-1])
        else:
            out.append(piece)
    return ''.join(out)


def extract_info(amplicon):
    name = amplicon.seq_id
    length = amplicon.length
    start = 0
    stop = 0
    for primer in amplicon.binding_sites:
        index, mismatches = amplicon.binding_sites[primer]
        if index > start:
            start, mis_start, primer_start = index, mismatches, primer
        if index < stop:
            stop, mis_stop, primer_stop = index, mismatches, primer
    primer_start = psearch_primer_to_iupac(primer_start)
    primer_stop = psearch_primer_to_iupac(primer_stop)
    return (name,
            start, stop,
            mis_start, mis_stop,
            primer_start, primer_stop)

def mung_file(handle):
    records = []
    amplifiers = read(handle).amplifiers
    for primer_set in amplifiers:
        for amplicon in amplifiers[primer_set]:
            records.append([primer_set] + list(extract_info(amplicon)))
    return DataFrame(records, columns=('primer_set', 'template', 'start', 'stop',
                                       'mismatch_start', 'mismatch_stop',
                                       'primer_start', 'primer_stop'))

def main():
    mung_file(sys.stdin).to_csv(sys.stdout, sep='\t', index=False)

if __name__ == '__main__':
    main()
