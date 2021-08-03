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
    description="takes read and index files from Illumina paired-end seq,  \
    demultiplexes the reads based on index, separates unknown index or low quality  \
    reads, and separates index-hopped reads"
)
parser.add_argument("-r1", help="read 1 file (fastq.gz)")
parser.add_argument("-r2", help="read 2 file (fastq.gz)")
parser.add_argument("-i1", help="index 1 file (fastq.gz)")
parser.add_argument("-i2", help="index 2 file (fastq.gz)")
parser.add_argument("-b", help="barcodes file (tsv)")
parser.add_argument("-l", help="length of each read")
parser.add_argument("-n", help="number of reads in the file")
args = parser.parse_args()

# define globals
read1_file: str = str(args.r1)
read2_file: str = str(args.r2)
index1_file: str = str(args.i1)
index2_file: str = str(args.i2)
barcodes_file: str = str(args.b)
read_length: int = int(args.l)
num_reads: int = int(args.n)

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

# dynamic file opening
####################################################################################################
# TODO: need to figure out how to open each of these files;
    # maybe we can have a list of unique variables that we index
    # in order to name each of the files, for example, if we
    # have a list of random character variable names and 26 input files
    # then we would index 0-25 of our random character variable list
barcodes = {}
with open(barcodes_file) as b:
    for line in b:
        line = line.rstrip().split('\t')
        barcodes[line[1]] = [line[0]]
# barcodes = sorted(barcodes)
print(barcodes)


# generate and open output files
# TODO: THIS DYNAMICALLY CREATES ALL OF OUR BARCODE FILES
    # SO CAN WE NOW JUST USE THE ITEMS IN THE LIST TO ACCESS
    # AND WRITE TO EACH OF THOSE FILES?
# TODO: SWITCH MY BARCODES DICTIONARY SO WE CAN USE THE 
    # BARCODE TO ACCESS THE CORRECT POSITION OF OUR LIST

# so here we are dynamically creating our output text files based on the names
    # of the barcodes that come in. we store the object in the dictionary
    # which can be accessed using the BARCODE. <- by creating access using
    # the BARCODE, then we can easily compare to our INDEXES that we find up above
        
        # FOR EXAMPLE:
            # if the index is a match, then write to the file that is associated
                # with that index

for f in barcodes:
    file = '{}.txt'.format(barcodes[f][0])
    barcodes[f].append(open(file, 'w'))
print(barcodes)

##################################################################################################

# open standard input and output files
with (gzip.open(read1_file) as rf1,
    gzip.open(read2_file) as rf2,
    gzip.open(index1_file) as if1,
    gzip.open(index2_file) as if2,
    open('read1_lowqual.fastq','w') as lowqual1,
    open('read2_lowqual.fastq','w') as lowqual2,
    open('read1_mismatch.fastq','w') as mismatch1,
    open('read2_mismatch.fastq','w') as mismatch2):

    # read 1st records from file
    rf1_record = [rf1.readline().decode('utf-8').rstrip() for _ in range(4)]
    rf2_record = [rf2.readline().decode('utf-8').rstrip() for _ in range(4)]
    if1_record = [if1.readline().decode('utf-8').rstrip() for _ in range(4)]
    if2_record = [if2.readline().decode('utf-8').rstrip() for _ in range(4)]

    # recordcount tracker
    recordcount = 1

    # loop and break at end of file
    while rf1_record[0] != '':

        # reverse complement index 2
        i2_rev = rev_complement(if2_record[1])

        # check quality score reads
        # TODO: may be a good idea to work a quality calculation check in here
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

        # check for unknown
        if 'N' in if1_record[1] or 'N' in if2_record[1] or low_qscore < 20:
            print('UNKNOWN OR LOW QUALITY!', if1_record[1], i2_rev)
            for component in rf1_record:
                lowqual1.write(component + '\n')
            for component in rf2_record:
                lowqual2.write(component + '\n')
        
        # # check for low quality
        # elif low_qscore < 20:
        #     print(new_header)
        #     print('LOW QUALITY!', if1_record[1], i2_rev)

        # check for index match
        elif if1_record[1] == i2_rev:
            print('MATCH!', if1_record[1], i2_rev)
            print('writing to {}.txt'.format(barcodes[if1_record[1]][0]))
            for component in if1_record:
                barcodes[if1_record[1]][1].write(component + '\n')
        else:
            print('MISMATCH!', if1_record[1], i2_rev)
            for component in rf1_record:
                mismatch1.write(component + '\n')
            for component in rf2_record:
                mismatch2.write(component + '\n')
        

        # read next record
        rf1_record = [rf1.readline().decode('utf-8').rstrip() for _ in range(4)]
        rf2_record = [rf2.readline().decode('utf-8').rstrip() for _ in range(4)]
        if1_record = [if1.readline().decode('utf-8').rstrip() for _ in range(4)]
        if2_record = [if2.readline().decode('utf-8').rstrip() for _ in range(4)]

        # increment recordcount and update statement
        recordcount += 1
        if recordcount % 10000000 == 0:
            print('Records processed:', recordcount)

# close each of our dynamic output files
for f in barcodes:
    barcodes[f][1].close()
    barcodes[f].pop(-1)

print(barcodes)




# Open each of the files to read through each corresponding record simultaneously.

#     For loop to read each record and compare indexes:
        
#         Reverse complement index2.
        
#         Calculate quality score for each of the current indexes.

#         Append indexes to header lines of read 1 and read 2.

#         If either index contains 'N':
#             store read record 1 to read1_lowqual category of dict
#             store read record 2 to read2_lowqual category of dict
#         elif either index is below quality threshold:
#             store read record 1 to read1_lowqual category of dict
#             store read record 2 to read2_lowqual category of dict
#         elif index1 != index2:
#             store read record 1 to read1_mismatch category of dict
#             store read record 2 to read2_mismatch category of dict
#         elif match:
#             store read record 1 to read1 for specific index category of dict
#             store read record 2 to read2 for specific index category of dict

# Output all files to specified output category (defines the file they end up in).
