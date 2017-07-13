#!/usr/bin/env python
"""
This script takes a fasta file and:

1) Removes any ambiguity codes (replaced by N)

2) Checks if the proportion of Ns is > user provided threshold

"""

import sys
import csv
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

def clean_fasta(f,o,c):
    seqs = SeqIO.parse(f,"fasta")
    sequences = []
    i = 1
    out_handle = open(o,'w')
    for seq in seqs:
        sequence = str(seq.seq)
        sequence = sequence.upper()
        # Recode anything that isn't A,C,T,G or N as N
        sequence = re.sub('[^ACTGN]','N',sequence)
        # Check coverage.
        N=sequence.count("N")
        Nperc=(float(N)/len(sequence))*100
        #print N,len(sequence),float(N)/len(sequence)
        if Nperc > int(c):
            sys.stderr.write("File: "+f+" -- "+
                             str(round(Nperc,2))+"% Ns is more than "+str(c)+\
                             "% Ns. Exiting\n")
            # Removes the empty file.
            os.remove(o)
            sys.exit()
        else:
            seq.seq = Seq(sequence,alphabet=generic_dna)
            SeqIO.write(seq,out_handle,"fasta")
    out_handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-f',
                       help = 'fasta file')

    parser.add_argument('-o',
                       help = 'output fasta file')
    
    parser.add_argument('-c',default = 20,
                       help = 'Completness threshhold (default 100%). '+\
                        'If more than c%'+\
                        'of the sequence is Ns then quit and give warning')

    args = vars(parser.parse_args())

    try:
        clean_fasta(args['f'],args['o'],args['c'])

    except IOError:
        #stdout is closed, no point in continuing
        #close explicitly to prevent problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass
