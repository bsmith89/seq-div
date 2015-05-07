#!/usr/bin/env python3
"Parse primersearch output files into tables."

import sys
from Bio.Emboss.PrimerSearch import read
from pandas import DataFrame

def parse_rec(hit_info):
    words = hit_info.split()
    name = words[0]
    start = int(words[6])
    stop = 0 - int(words[15][1:-1])
    mis_start = int(words[8])
    mis_stop = int(words[17])
    primer_start = words[1]
    primer_stop = words[10]
    return (name, start, stop, mis_start, mis_stop, primer_start, primer_stop)

def mung_file(handle):
    records = []
    amplifiers = read(handle).amplifiers
    for primer_set in amplifiers:
        for amplicon in amplifiers[primer_set]:
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
            records.append((name, primer_set,
                            start, stop,
                            mis_start, mis_stop,
                            primer_start, primer_stop))
    return DataFrame(records, columns=('template', 'primer_set', 'start', 'stop',
                                       'mismatch_start', 'mismatch_stop',
                                       'primer_start', 'primer_stop'))

def main():
    mung_file(sys.stdin).to_csv(sys.stdout, sep='\t', index=False)

if __name__ == '__main__':
    main()
