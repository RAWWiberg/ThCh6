#!/usr/bin/env python
"""
This script takes a multi-sequence fasta file and user suplied list of IDs and
will print only those sequences that have the supplied IDs

Fasta headers should have the format: species_id|gene_id
"""

import csv
import argparse
import sys
import os
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord

def read_list(fil):
    fil_reader = csv.reader(fil,delimiter = "\t")
    IDlist = [line[0] for line in fil_reader]
    return IDlist

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-IDs',
                       help = 'single ID to extract or list of IDs in a file')

    parser.add_argument('-splitID',
                       help = 'split (1) or not (0) the sequence IDs')

    parser.add_argument('-fasta',
                       help = 'fasta filename')

    parser.add_argument('-task',default = 'f',
                       help = 'f = (def) filter the fasta to include only IDs,'+\
                        'fs = filter and split fasta into one file per ID')

    parser.add_argument('-outfile',
                       help = 'output filename, only if task = f')

    parser.add_argument('-outdir',
                       help = 'output directory')

    args = vars(parser.parse_args())

    try:
        try:
            fil = open(args['IDs'],'r')
            IDs = read_list(fil)
            #are there any duplicates in the list?
            if len(IDs) != len(set(IDs)):
                print "There are duplicates in the list: \n"
            #print "\nThese are the IDs in the ID list file: ", IDs, "\n\n"# script tester line
            print "\nSequences in IDs file: ",len(IDs),"\n"
            fastafile = SeqIO.parse(open(args['fasta'],'r'),"fasta")
            i = 0
            outdir = args['outdir']
            if args['splitID'] == '1':
                if args['task'] == 'fs':
                    for sequence in fastafile:
                        seqid = sequence.id.split("|")[1]
                        if any(ID == seqid for ID in IDs):
                            ID = [iD for iD in IDs if iD == seqid]
                            output_handle = open(outdir+"/"+"".join(ID)+".fasta",'w')
                            SeqIO.write(sequence,output_handle,"fasta")
                            output_handle.close()
                            #print sequence.format("fasta")
                            i = i+1        
                elif args['task'] == 'f':
                    output_handle = open(outdir+"/"+args['outfile'],'w')
                    for sequence in fastafile:
                        seqid = sequence.id.split("|")[1]
                        if any(ID == seqid for ID in IDs):
                            ID = [iD for iD in IDs if iD == seqid]
                            SeqIO.write(sequence,output_handle,"fasta")
                            #print sequence.format("fasta")
                            i = i+1
                    output_handle.close()
            else:
                if args['task'] == 'fs':
                    for sequence in fastafile:
                        seqid = sequence.id
                        if any(ID == seqid for ID in IDs):
                            ID = [iD for iD in IDs if iD == seqid]
                            output_handle = open(outdir+"/"+"".join(ID)+".fasta",'w')
                            SeqIO.write(sequence,output_handle,"fasta")
                            output_handle.close()
                            #print sequence.format("fasta")
                            i = i+1        
                elif args['task'] == 'f':
                    output_handle = open(outdir+"/"+args['outfile'],'w')
                    for sequence in fastafile:
                        seqid = sequence.id
                        if any(ID == seqid for ID in IDs):
                            ID = [iD for iD in IDs if iD == seqid]
                            SeqIO.write(sequence,output_handle,"fasta")
                            #print sequence.format("fasta")
                            i = i+1
                    output_handle.close()
                    
            print "\nSequences processed: ",i,"\n"

        except IOError:
            #The input here is not a filename, do something else
            #print "This is not a filename, must be a single ID"#script tester line
            ID = args['IDs']
            fastafile = SeqIO.parse(open(args['fasta'],'r'),"fasta")
            outdir = args['outdir']
            if args['splitID'] == "1":
                for sequence in fastafile:
                    seqid = sequence.id.split("|")[1]
                    if ID == seqid:
                        output_handle = open(outdir+"/"+args['outfile'],'w')
                        SeqIO.write(sequence,output_handle,"fasta")
                        output_handle.close()
                        #print sequence.format("fasta")
            else:
                for sequence in fastafile:
                    seqid = sequence.id
                    if ID == seqid:
                        output_handle = open(outdir+"/"+args['outfile'],'w')
                        SeqIO.write(sequence,output_handle,"fasta")
                        output_handle.close()
                        #print sequence.format("fasta")                

    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass
