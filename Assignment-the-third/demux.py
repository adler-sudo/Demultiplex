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
parser.add_argument("-c", help="index mean quality cutoff")
args = parser.parse_args()

# define globals
read1_file: str = str(args.r1)
read2_file: str = str(args.r2)
index1_file: str = str(args.i1)
index2_file: str = str(args.i2)
barcodes_file: str = str(args.b)
qual_cutoff: int = int(args.c) # TODO: put this in parser

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

# establish barcodes dictionary of dictionary (key = index, value = dict(barcode, read1file, read2file, counter))
barcodes = {}
with open(barcodes_file) as b:
    for line in b:
        line = line.rstrip().split('\t')
        
        # generate comprehensive 
        barcodes[line[1]] = {
            'barcode':line[0],
            'read1file':'blank',
            'read2file':'blank',
            'counter':0
        }

# so here we are dynamically creating our output text files based on the names
for f in barcodes:
    read1_match = '{}_read1.fastq'.format(barcodes[f]['barcode'])
    read2_match = '{}_read2.fastq'.format(barcodes[f]['barcode'])
    barcodes[f]['read1file'] = open(read1_match, 'w')
    barcodes[f]['read2file'] = open(read2_match, 'w')

# generate swap dictionary and counter for tracking instance of each swap
swap_counter = 0
swap_dict = {}
for i1 in barcodes:
    for i2 in barcodes:
        if i1 != i2:
            index_combo = i1 + '-' + i2
            swap_dict[index_combo] = 0
print('length of swap dict:',len(swap_dict))

# open standard input files
rf1 = gzip.open(read1_file, 'rt')
rf2 = gzip.open(read2_file, 'rt')
if1 = gzip.open(index1_file, 'rt')
if2 = gzip.open(index2_file, 'rt')

# open standard output files
lowqual1 = open('read1_lowqual.fastq','w')
lowqual2 = open('read2_lowqual.fastq','w')
mismatch1 = open('read1_mismatch.fastq','w')
mismatch2 = open('read2_mismatch.fastq','w')

# recordcount tracker
recordcount = 0
unknown_counter = 0
lowqual_counter = 0
match_counter = 0

# loop and break at end of file
while True:

    # read next record
    rf1_record = [rf1.readline().rstrip() for _ in range(4)]
    rf2_record = [rf2.readline().rstrip() for _ in range(4)]
    if1_record = [if1.readline().rstrip() for _ in range(4)]
    if2_record = [if2.readline().rstrip() for _ in range(4)]

    # break at end of file
    if rf1_record[0] == '':
        break

    # increment recordcount and update statement
    recordcount += 1
    if recordcount % 100000 == 0:
        print('Records processed:', recordcount)

    # reverse complement index 2
    i2_rev = rev_complement(if2_record[1])

    # generate new header
    index_combo = if1_record[1] + '-' + i2_rev
    new_header1 = rf1_record[0] + ':' + index_combo
    new_header2 = rf2_record[0] + ':' + index_combo
    rf1_record[0] = new_header1 # TODO: probably a more creative way to assign the new header?
    rf2_record[0] = new_header2

    # if contains unknown base in index, write to low qual
    if 'N' in if1_record[1] or 'N' in if2_record[1]:
        for component in rf1_record:
            lowqual1.write(component + '\n')
        for component in rf2_record:
            lowqual2.write(component + '\n')
        unknown_counter += 1
        continue

    # if not valid barcode, write to low qual
    if if1_record[1] not in barcodes or i2_rev not in barcodes:
        for component in rf1_record:
            lowqual1.write(component + '\n')
        for component in rf2_record:
            lowqual2.write(component + '\n')
        lowqual_counter += 1
        continue

    # check quality score reads
    if1_qscore = sum([int(convert_phred(letter)) for letter in if1_record[3]]) / len(if1_record[3])
    if2_qscore = sum([int(convert_phred(letter)) for letter in rf2_record[3]]) / len(rf2_record[3])

    # lowest qscore
    low_qscore = min([if1_qscore,if2_qscore])

    # check for N, unknown barcode, or low quality score and write to low_qual if so
    if low_qscore < qual_cutoff:
        for component in rf1_record:
            lowqual1.write(component + '\n')
        for component in rf2_record:
            lowqual2.write(component + '\n')

    # check for index match and write to match files for corresponding index
    elif if1_record[1] == i2_rev:
        for component in rf1_record:
            barcodes[if1_record[1]]['read1file'].write(component + '\n')
        for component in rf2_record:
            barcodes[if1_record[1]]['read2file'].write(component + '\n') # added reference to second read file (added post 8/2 run)
        barcodes[if1_record[1]]['counter'] += 1
        match_counter += 1

    # otherwise write to mismatch file
    else:
        for component in rf1_record:
            mismatch1.write(component + '\n')
        for component in rf2_record:
            mismatch2.write(component + '\n')
        
        # increment swap dict
        swap_dict[index_combo] += 1
        
        # increment swap counter
        swap_counter += 1

# close gzip files
rf1.close()
rf2.close()
if1.close()
if2.close()

# close standnard output files
lowqual1.close()
lowqual2.close()
mismatch1.close()
mismatch2.close()

# close each of our dynamic output files
for f in barcodes:
    barcodes[f]['read1file'].close()
    barcodes[f]['read2file'].close() # added post 8/2 run

# summary barcodes
sum_barcodes = {}

# create summary file
with open('summary_stats.md','w') as s:
    
    # calculate overall stats
    match_perc_total = round(match_counter / recordcount * 100, 2)
    swap_perc_total = round(swap_counter / recordcount * 100, 2)
    unknown_perc_total = round(unknown_counter / recordcount * 100, 2)
    lowqual_perc_total = round(lowqual_counter / recordcount * 100, 2)

    # calculate percentage of each sample
    for f in barcodes:
        barcode_id = barcodes[f]['barcode']
        perc_reads = round(barcodes[f]['counter'] / recordcount * 100, 2)
        sum_barcodes[barcode_id] = perc_reads
    
    # write description
    s.write('# Summary Stats\n')
    s.write('## Overall\n')
    s.write('Number of records processed: ' + str(recordcount) + '<br>')
    s.write('Index mean quality score cutoff: ' + str(qual_cutoff) + '<br>')
    s.write('Percentage correctly indexed: ' + str(match_perc_total) + '%' + '<br>')
    s.write('Percentage of index swaps: ' + str(swap_perc_total) + '%' + '<br>')
    s.write('Percentage unknown index: ' + str(unknown_perc_total) + '%' + '<br>')
    s.write('Percentage low quality index: ' + str(lowqual_perc_total) + '%' + '<br>')

    # sort and write percentage and adapter
    s.write('\n## Percentage of records from each adapter:\n')
    s.write('```\n')
    s.write('adapter' + '\t' + 'perc_records\n')
    for sb in sorted(sum_barcodes):
        s.write('{}'.format(sb) + '\t' + str(sum_barcodes[sb]) + "%" + "\n")
    s.write('```')

    # write percentage for each index swap
    s.write('\n## Reads with each potential index swap:\n')
    s.write('```\n')
    s.write('index_combo' + '\t' + 'count' + '\t' + 'perc_of_all_mismatches\n')
    for swap in swap_dict:
        swap_perc = round(swap_dict[swap] / swap_counter * 100, 2)
        s.write(swap + '\t' + str(swap_dict[swap]) + '\t' + str(swap_perc) + '%' + '\n')
    s.write('```')

    
