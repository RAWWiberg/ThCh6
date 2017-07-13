#!/usr/bin/env python
"""
This script takes a phylip file and checks that the formatting is
appropriate for paml:

Checks sequence name length.
Checks length divisible by 3 (i.e. there is a full ORF)

if all above checks pass:
Also removes terminal STOP codons. 
"""

import sys
import argparse
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Gapped, generic_dna
import re


def pamlcleaner(phy):
    outputfile = open(phy.replace('.phy', '_clean.phy'), 'w')
    alignments = open(phy, 'rb')
    content = alignments.readlines()
    # reads in contents of phylip file to a list, if file is too large
    # this might be a problem
    content = [i.replace('\n', '') for i in content]
    d = []
    d2 = []
    head = content[0].split(' ')
    seqs = content[1:] # everything in file except header
    for line in seqs:
        #collect names of sequences
        if len(line) < 30 and any(i.islower() for i in line):
            # a name line in the phylip must have less than 30
            # characters and contain at least one lower case character
            # sequence lines must not fit this description for this to work
            d.append(line)       
    seqs = re.split('|'.join(d), ''.join(seqs))
    # splits the sequences into single strings based on names in d
    # output is list of sequences in same order as names in d
    seqs = [Seq(seq, alphabet = Gapped(generic_dna, "-")) for seq in seqs[1:]]
    # convert each string into a seq object which can be checked by
    # BioPython functions
    #print len(seqs), int(head[0]) # script tester line
    if len(seqs) != int(head[0]):
        # number of sequences matches header
        sys.stderr.write("more sequences than mentioned in file header: "+
                         phy+" EXITING!\n")
        sys.exit()

    name = 0
    for seq in seqs:
        # I should let PAML warn about stop codons as this ungapping
        # step produces weird sequences
        #if str(seq.ungap().translate(1))[:-1].find("*") > 0:
        #    # premature stop codons found
        #    sys.stderr.write("premature stop codon found in: "+phy+
        #                    " - "+d[name]+" EXITING!\n")
        #    sys.exit()
        if len(seq) != int(head[1]) or (int(head[1])-3)%3 != 0:
            # sequence not same length as described in header
            sys.stderr.write("sequence lengths do not correspond"+
                                " to header OR not a multiple of 3 "+
                                ": "+phy +" - "+d[name]+" EXITING!\n")
            sys.exit()
        name = name + 1
    # all checks have passed, print the sequences (minus
    # last three bases (terminal STOP codon))
    length = int(head[1])
    head[1] = str(int(head[1])-3)
    outputfile.write(' '.join(head)+'\n')
    for i in range(0,len(d)):
        outputfile.write(d[i]+'\n')
        lines = length/60
        s = 0
        e = 60
        k = 1
        while k <= lines:
            outputfile.write(str(seqs[i])[s:e]+'\n')
            k = k + 1
            s = e
            e = e + 60
        outputfile.write(str(seqs[i])[s:-3]+'\n')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-phy',
                       help = 'phylip format file for checking')
    
    args = vars(parser.parse_args())

    try:
        phy = args['phy']
        pamlcleaner(phy)
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass


