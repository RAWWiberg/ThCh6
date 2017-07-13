#!/usr/bin/python
"""
This script takes the directory fasta files
and strips name from each sequence entry:

e.g.:
>sp_name|gene_model_1 --> >sp_name

Assumes every fasta in the directory is formatted the same way.
"""

import sys
import csv
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def strip_names(fas,out_handle):
    outfil = open(out_handle,'w')
    seqs = SeqIO.parse(fas,"fasta")
    sequences = []
    i = 1
    for seq in seqs:
        #print seq.id, len(seq.id.split("|")) > 1
        if len(seq.id.split("|")) > 1:
            seq.id = seq.id.split("|")[0]
            seq.name = ""
            seq.description = ""
            sequences.append(seq)
        else:
            seq.id = "sp_" + str(i)
            seq.name = ""
            seq.description = ""
            sequences.append(seq)
            i = i + 1
    SeqIO.write(sequences, outfil, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-f',
                       help = 'input fasta file')

    parser.add_argument('-of',
                       help = 'output file handle')

    args = vars(parser.parse_args())

    try:
        fas = open(args['f'],'r')
        out_handle = args['of']+'.fasta'
        strip_names(fas,out_handle)
                
    except IOError:
        #stdout is closed, no point in continuing
        #close explicitly to prevent problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass
