#!/usr/bin/env python
"""
This script takes a multi-sequence fasta file filters sequences on:

length
completeness

Tested for:
ensembl format
gffread (-x and -y options) format

Can also try to print only the longest CDS where several exist for the same gene
if there are ties, one will be picked at random.
"""

import csv
import argparse
import sys
import textwrap
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord


def fa_nt_filter(fa_fil, min_length, compl,lon, seqtype):
    fas = SeqIO.parse(fa_fil, "fasta")
    min_length=int(min_length)
    s = 1
    if lon == '1':
        #try to print only longest CDS
        #need to make a dictionary of genes and transcript lengths.
        seq_dict = {}
        if seqtype == 'nucl':
            for sequence in fas:
                name = sequence.id
                descr = sequence.description
                descr = descr.split(" ")
                print [i for i in descr if 'gene:' in i]
                print "SEQUENCE: ",s,len(sequence.seq),"\n",name,"\n",descr
                s = s + 1
        elif seqtype == 'pep':
            for sequence in fas:
                name = sequence.id
                descr = sequence.description
                descr = descr.split(" ")
                #look for ensembl style annotation in headers
                gene_id = "".join([i for i in descr if 'gene:' in i])
                trans_id = "".join([i for i in descr if 'transcript:' in i])
                gene_id = gene_id.split(":")[1]
                trans_id = trans_id.split(":")[1]
                if len(sequence.seq) >= length:
                    print "SEQUENCE: ",s,len(sequence.seq),"\n\n",
                    print gene_id, trans_id
                    print name,"\n\n",descr
                    s = s + 1
    elif lon == '0':
        if seqtype == 'nucl':
            for sequence in fas:
                name = sequence.id
                descr = sequence.description
                descr = descr.split(" ")
                if len(sequence.seq) > min_length:
                    print ">"+sequence.name
                    pr_seq=textwrap.wrap(str(sequence.seq),width=60)
                    print "\n".join(pr_seq)
                    print "\n"
                    #print "SEQUENCE: ",s,len(sequence.seq),"\n",name,"\n",descr
                    s = s + 1
        elif seqtype == 'pep':
            for sequence in fas:
                name = sequence.id
                descr = sequence.description
                descr = descr.split(" ")
                if len(sequence.seq) > min_length:
                    print ">"+sequence.name
                    pr_seq=textwrap.wrap(str(sequence.seq),width=60)
                    print "\n".join(pr_seq)
                    print "\n"
                    #print "SEQUENCE: ",s,len(sequence.seq),"\n",name,"\n",descr
                    s = s + 1

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-f',
                        help = 'input fasta file.')

    parser.add_argument('-min_length', default = 1,
                        help = 'minimum length of sequence (nucl/amino acids).')

    parser.add_argument('-min_compl',default = 0.8,
                        help = 'minimum completeness, proportion of sites'+\
                        'that are allowed to be missing.')

    parser.add_argument('-type', default = 'nucl',
                        help = 'sequence type, nucleotide (nucl) or peptide'+\
                        'pep.')

    parser.add_argument('-long', default = '0',
                        help = 'if (1) will try to output only longest of multiple'+\
                        'CDSs, if (0) will output everything.')

    args = vars(parser.parse_args())


    try:
        fil = open(args['f'],'r')
        fa_nt_filter(fil,args['min_length'],args['min_compl'],
                     args['long'],args['type'])

    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass

