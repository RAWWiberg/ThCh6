#!/usr/bin/env python
"""
This script takes a list of fasta files and concatenates the sequences for each
entry.

Entries must have the same names.
"""

import csv
import argparse
import sys
from glob import glob
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def concSeqs(f,seq_dict):
    fas = SeqIO.parse(f, "fasta")
    for record in fas:
        if record.id in seq_dict.keys():
            # The id already exists in the dictionary
            seq_dict[record.id] = seq_dict[record.id]+str(record.seq)
        else:
            # The id doesn't exist
            seq_dict[record.id]=str(record.seq)
            
    return seq_dict

def printDict(seq_dict):
    for record in seq_dict.keys():
        seq_record = SeqRecord(Seq(seq_dict[record]),id=record,
                               description="")
        print seq_record.format('fasta')

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])
    
    parser.add_argument('-f',nargs="+",
                        help = 'fasta file')

#    parser.add_argument('-o',
 #                       help = 'out file')
#
    args = vars(parser.parse_args())


# Start a dictionary of seq objects
seq_dict={}
# For each sequence file in the list of files
for f in args['f']:
    #print f
    concSeqs(f,seq_dict)
printDict(seq_dict)


