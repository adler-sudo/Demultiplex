#!/usr/bin/env python

# import modules
import gzip
import Bioinfo
import argparse
import os
import gzip
import numpy as np
import matplotlib.pyplot as plt
import io
from Bioinfo import convert_phred


# define arparser
parser = argparse.ArgumentParser(
    description="takes barcodes, read files, and index files from Illumina paired-end  \
    seq, demultiplexes the reads based on index, separates unknown index or low quality  \
    reads, and separates index-hopped reads"
)
parser.add_argument("-r1", help="read 1 file (fastq.gz)")
parser.add_argument("-r2", help="read 2 file (fastq.gz)")
parser.add_argument("-i1", help="index 1 file (fastq.gz)")
parser.add_argument("-i2", help="index 2 file (fastq.gz)")
parser.add_argument("-b", help="barcodes file (tsv)")
args = parser.parse_args()

# define globals
read1_file: str = str(args.r1)
read2_file: str = str(args.r2)
index1_file: str = str(args.i1)
index2_file: str = str(args.i2)
barcodes_file: str = str(args.b)

# file format check
assert read1_file.endswith('.fastq.gz'), "Read 1 file must be .fastq.gz"
assert read2_file.endswith('.fastq.gz'), "Read 2 file must be .fastq.gz"
assert index1_file.endswith('.fastq.gz'), "Index 1 file must be .fastq.gz"
assert index2_file.endswith('.fastq.gz'), "Index 2 file must be .fastq.gz"
assert barcodes_file.endswith('.tsv'), "Barcodes file must be .tsv"

# define reverse function
def rev_complement(strand):
    """returns the reverse complement strand of the string"""

    # reverse strand
    rev_strand = strand[::-1]

    # create the dictionary for variable assignment
    a = {
        'G':'C',
        'C':'G',
        'A':'T',
        'T':'A',
        'N':'N' # this will work for now but not sure if this is agnostic
    }
    
    # create and return the compliment strand
    rev_comp = ''
    for i in rev_strand:
        complement = a[i]
        rev_comp = rev_comp + complement
    
    return rev_comp

# dynamically open a read1 and read2 for each barcode
barcodes = {}
with open(barcodes_file) as b:
    for line in b:
        line = line.rstrip().split('\t')
        barcodes[line[1]] = [line[0]]

# so here we are dynamically creating our output text files based on the names
for f in barcodes:
    read1_match = '{}_read1.fastq.gz'.format(barcodes[f][0])
    read2_match = '{}_read2.fastq.gz'.format(barcodes[f][0]) # added this post 8/2 run
    barcodes[f].append(gzip.open(read1_match, 'wb'))
    barcodes[f].append(gzip.open(read2_match, 'wb')) # added this post 8/2 run

# open standard input and output files
with gzip.open(read1_file, 'rt') as rf1, \
    gzip.open(read2_file, 'rt') as rf2, \
    gzip.open(index1_file, 'rt') as if1, \
    gzip.open(index2_file, 'rt') as if2, \
    gzip.open('read1_lowqual.fastq.gz','wb') as lowqual1, \
    gzip.open('read2_lowqual.fastq.gz','wb') as lowqual2, \
    gzip.open('read1_mismatch.fastq.gz','wb') as mismatch1, \
    gzip.open('read2_mismatch.fastq.gz','wb') as mismatch2:

    # read 1st records from file
    rf1_record = [rf1.readline().rstrip() for _ in range(4)]
    rf2_record = [rf2.readline().rstrip() for _ in range(4)]
    if1_record = [if1.readline().rstrip() for _ in range(4)]
    if2_record = [if2.readline().rstrip() for _ in range(4)]

    # recordcount tracker
    recordcount = 1

    # loop and break at end of file
    while rf1_record[0] != '':

        # reverse complement index 2
        i2_rev = rev_complement(if2_record[1])

        # check quality score reads
        # TODO: may be a good idea to work a quality calculation check in here (unittest)
        # TODO: turn this into a function
        rf1_qscore = sum([int(convert_phred(letter)) for letter in rf1_record[3]]) / len(rf1_record[3])
        rf2_qscore = sum([int(convert_phred(letter)) for letter in rf2_record[3]]) / len(rf2_record[3])
        if1_qscore = sum([int(convert_phred(letter)) for letter in if1_record[3]]) / len(if1_record[3])
        if2_qscore = sum([int(convert_phred(letter)) for letter in rf2_record[3]]) / len(rf2_record[3])

        # lowest qscore
        low_qscore = min([rf1_qscore,rf2_qscore,if1_qscore,if2_qscore])

        # generate new header
        new_header1 = rf1_record[0] + ':' + if1_record[1] + '-' + i2_rev
        new_header2 = rf2_record[0] + ':' + if1_record[1] + '-' + i2_rev
        rf1_record[0] = new_header1 # TODO: probably a more creative way to assign the new header?
        rf2_record[0] = new_header2

        # check for N, unknown barcode, or low quality score and write to low_qual if so
        if 'N' in if1_record[1] or 'N' in if2_record[1] or if1_record[1] not in barcodes or low_qscore < 20:
            for component in rf1_record:
                lowqual1.write((component + '\n').encode())
            for component in rf2_record:
                lowqual2.write((component + '\n').encode())

        # check for index match and write to match files for corresponding index
        elif if1_record[1] == i2_rev:
            for component in rf1_record:
                barcodes[if1_record[1]][1].write((component + '\n').encode())
            for component in rf2_record:
                barcodes[if1_record[1]][2].write((component + '\n').encode()) # added reference to second read file (added post 8/2 run)
        
        # otherwise write to mismatch file
        else:
            for component in rf1_record:
                mismatch1.write((component + '\n').encode())
            for component in rf2_record:
                mismatch2.write((component + '\n').encode())
        
        # read next record
        rf1_record = [rf1.readline().rstrip() for _ in range(4)]
        rf2_record = [rf2.readline().rstrip() for _ in range(4)]
        if1_record = [if1.readline().rstrip() for _ in range(4)]
        if2_record = [if2.readline().rstrip() for _ in range(4)]

        # increment recordcount and update statement
        recordcount += 1
        if recordcount % 1000000 == 0:
            print('Records processed:', recordcount)

# close each of our dynamic output files
for f in barcodes:
    barcodes[f][1].close()
    barcodes[f][2].close() # added post 8/2 run
    
    # get rid of file instances following close
    barcodes[f].pop(-1)
    barcodes[f].pop(-1)
