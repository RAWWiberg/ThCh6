#!/usr/bin/env python
"""
This script takes a fasta file and re-names the sequences
in a simple integer order with a common attached handle provided by
user (default).

e.g. user provides handle "dmel"

script changes all individual sequences to
dmel|dmel000001
dmel|dmel000002
dmel|dmel000003
...

Alternative 1: The script can simply add the handle to the beginning of the
existing sequence ID

Alternative 2: The script changes the sequence ID to the file ID (up to the
first "_" or space " " and adds the handle to the beginning

Also produces a file called [handle]_names_file.txt with a table of new
names and the old matching names.
"""

import csv
import argparse
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord

def fa_modifier(fa_fil, handle, a):
    fas = SeqIO.parse(fa_fil, "fasta")
    i = 1
    if a == '0':
        outfil = csv.writer(open(handle+"_names_file.txt",'w'), delimiter = "\t")
        for sequence in fas:
            old_nam = sequence.id
            new_nam = handle + "|" + handle + str(i).zfill(5)
            sequence.id = new_nam
            sequence.name = new_nam
            sequence.description = ""
            i = i + 1
            outfil.writerow([old_nam, new_nam])
            print sequence.format('fasta')
    elif a == '1':
        for sequence in fas:
            old_nam = sequence.id
            new_nam = handle + "|" + sequence.id
            sequence.name = new_nam
            sequence.id = new_nam
            sequence.description = ""
            print sequence.format('fasta')
    elif a == '2':
        for sequence in fas:
            #print fa_fil.replace('_',' ').split(' ')[0]
            old_nam = sequence.id
            new_nam = handle + "|" + fa_fil.replace('_',' ').split(' ')[0]
            sequence.name = new_nam
            sequence.id = new_nam
            sequence.description = ""
            print sequence.format('fasta')

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-f',
                        help = 'fasta file')

    parser.add_argument('-ha',
                        help = 'user provided handle')

    parser.add_argument('-a', default = 0,
                        help = 'Alternative: change gene names to numers \
                        (0), just add handle (1), change to filename and add \
                        handle (2). Default = 0')

    args = vars(parser.parse_args())



fa_modifier(args['f'],args['ha'], args['a'])
