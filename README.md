# Thesis Chapter 6

This repository contains a description of the steps that were taking in the processing of genomic samples and subsequent 
analysis. It contains the following files.


# 1) Pipeline descriptions
These are descriptions of the bioinformatic steps taken in the processing of reads,
mapping, and analysis

#1.1) mapping_pipeline.txt

#1.2) BUSCO_to_phylo_pipeline.txt

#1.3) PAML_pipeline.txt

#1.4) CAFE_pipeline.txt

#1.5) ANGSD_pipeline.txt


# 2) R scripts
These are R scripts for the analysis and plotting of data

#2.1) crows_orthologs.R

#2.2) crows_BUSCO_results_analysis.R

#2.3) crows_paml_results_analysis.R

#2.4) crows_cwoo_cmon_thetasFST.R

#2.5) crows_pop_gen_sims.R

#2.6) sfs_stats.R

#2.7) crows_tree.R


# 3) Python scripts

These are python scripts to perform some bioinformatic processing and analysis.
They were all authored by me

#3.1) fastaStats.py

Simply counts the numbers of different nucleotides as well as Ns

#3.2) getFasta.py

Takes a list of sequence IDs and a fasta file with multiple sequences in and outputs
only those sequences in the list.

#3.3) cleanFasta.py

"Cleans" a fasta sequence by replacing any ambiguity codes with Ns. Will throw
a warning and exit if the proportion of Ns is greater than a user specified value

#3.3) pamlCleaner.py

"Cleans" a paml input file. Makes sure sequence is multiple of 3 in length, removes
terminal stop codons.

#3.3) createCodemlCtlFile.py

Creates a paml "control" file.

#3.4) codemlResultsParser.py

Parses a codeml output file to extract results into a single line.

#3.5) concatFasta.py

Concatenates the given fasta files

#3.6) renameFasta.py

Adds a new handle to the start of each sequence name

#3.7) stripFastaNames.py

Removes the ending of each header in a multi sequence fasta file to keep only the first
element of the header.

#3.8) filterFasta.py

Filters a multisequence fasta file



