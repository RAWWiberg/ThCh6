#!/usr/bin/env python
"""
This script takes a fasta file and for each entry outputs:
length (total number of sites)
number "Ns"
number "As"
number "Ts"
number "Gs"
number "Cs"
number lowercase
number uppercase
"""

import csv
import argparse
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord

def fa_nt_p_stats(fa_fil):
    fas = SeqIO.parse(fa_fil, "fasta")
    print "\t".join(["sequence_id","length","upper_count","lower_count",
                     "N_count","A_count","T_count","G_count","C_count"])
    for sequence in fas:
        name = sequence.id
        length = str(len(sequence.seq))
        n_lower = str(sum([int(b.islower()) for b in sequence.seq]))
        n_upper = str(sum([int(b.isupper()) for b in sequence.seq]))
        n_N = str(str(sequence.seq).count("N")+str(sequence.seq).count("n"))
        n_A = str(str(sequence.seq).count("A")+str(sequence.seq).count("a"))
        n_T = str(str(sequence.seq).count("T")+str(sequence.seq).count("t"))
        n_G = str(str(sequence.seq).count("G")+str(sequence.seq).count("g"))
        n_C = str(str(sequence.seq).count("C")+str(sequence.seq).count("c"))
        printline = [name,length,n_upper,n_lower,n_N,n_A,n_T,n_G,n_C]
        print "\t".join(printline)

def fa_nt_o_stats(fa_fil):
    fas = SeqIO.parse(fa_fil, "fasta")
    print "\t".join(["n_seqs","length","upper_count","lower_count",
                     "N_count","A_count","T_count","G_count","C_count"])
    seqs = 0
    length = 0
    n_lower = 0
    n_upper = 0
    n_N = 0
    n_A = 0
    n_T = 0
    n_G = 0
    n_C = 0
    for sequence in fas:
        seqs = seqs + 1
        length = length+len(sequence.seq)
        n_lower = n_lower+sum([int(b.islower()) for b in sequence.seq])
        n_upper = n_upper+sum([int(b.isupper()) for b in sequence.seq])
        n_N = n_N+str(sequence.seq).count("N")+str(sequence.seq).count("n")
        n_A = n_A+str(sequence.seq).count("A")+str(sequence.seq).count("a")
        n_T = n_T+str(sequence.seq).count("T")+str(sequence.seq).count("t")
        n_G = n_G+str(sequence.seq).count("G")+str(sequence.seq).count("g")
        n_C = n_C+str(sequence.seq).count("C")+str(sequence.seq).count("c")

    printline = [str(seqs),str(length),str(n_upper),str(n_lower),str(n_N),str(n_A),str(n_T),str(n_G),str(n_C)]
    print "\t".join(printline)

def fa_p_p_stats(fa_fil):
    fas = SeqIO.parse(fa_fil, "fasta")
    print "\t".join(["sequence_id","length","stop","stop_count"])
    for sequence in fas:
        name = sequence.id
        length = str(len(sequence.seq))
        stop = 0
        if '.' in sequence.seq:
            stop = 1
        n_stop = sequence.seq.count('.')
        printline = [name,length,str(stop),str(n_stop)]
        print "\t".join(printline)

def fa_p_o_stats(fa_fil):
    fas = SeqIO.parse(fa_fil, "fasta")
    print "\t".join(["n_seqs","length","stop","stop_count"])
    seqs = 0
    length = 0
    stop = 0
    n_stop = 0
    for sequence in fas:
        seqs = seqs + 1
        length = length+len(sequence.seq)
        if '.' in sequence.seq:
            stop = 1
        n_stop = n_stop + sequence.seq.count('.')
    printline = [str(seqs),str(length),str(stop),str(n_stop)]
    print "\t".join(printline)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-f',
                        help = 'input fasta file')
    
    parser.add_argument('-t', default = 'o',
                        help = 'o = overall counts, p = per sequence')

    parser.add_argument('-s', default = 'n',
                        help = 'n = nucleotide sequence, p = protein sequence')

    args = vars(parser.parse_args())

    try:
        if args['t'] == 'p':
            if args['s'] == 'n':
                fa_nt_p_stats(args['f'])
            if args['s'] == 'p':
                fa_p_p_stats(args['f'])
        elif args['t'] == 'o':
            if args['s'] == 'n':
                fa_nt_o_stats(args['f'])
            if args['s'] == 'p':
                fa_p_o_stats(args['f'])

        
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass

