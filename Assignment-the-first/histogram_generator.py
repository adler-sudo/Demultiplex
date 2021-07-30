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


# define arparser
parser = argparse.ArgumentParser(description="develop historgram of means for each position of a series of reads")
parser.add_argument("-f", help="input fastq file")
parser.add_argument("-l", help="length of each read")
parser.add_argument("-n", help="number of reads in the file")
parser.add_argument("-o", help="output histogram file name")
args = parser.parse_args()

# define globals
input_file: str = args.f
read_length: int = int(args.l)
num_reads: int = int(args.n)
output_file: str = args.o

# check if fastq file is gzipped
if input_file.endswith('.gz'):
    gzipped: bool = True
elif input_file.endswith('.fastq'):
    gzipped: bool = False
else:
    raise ValueError("input file must be .gz or .fastq file")

# define postion qscore average function
def qscore_average_function(file: str, read_length: int, num_reads: int, gzipped=True):
    """
    Given fatsq file or gzipped fastq file, generate an array of mean qscore by position.
    
    Paramteters:
    ------------
        file: fastq file to perfom calculations on
        
        read_length: length of each read in fastq file
        
        num_reads: number of reads in fastq file

        gzipped: is fastq file gzipped
    
    Returns:
    --------
        mean:
            - array of mean qscore where index is equal to position. 
            
        linecount:
            - total lines in file
    """

    # initiate score sum
    linecount = 0
    
    # initiate our array
    mean = np.zeros([read_length,1])
    
    # determine if gzipped
    if gzipped:

        with gzip.open(file, 'rb') as f:
        
            # loop through each line
            while True:
                newline = f.readline().rstrip()
                newline = newline.decode('utf-8')
                
                # break if empty (end of file)
                if newline == '':
                    break
                
                # update statement
                if linecount % 500000 == 0:
                    print(linecount)
                
                # grab quality score line of record
                if linecount % 4 == 3:
                    
                    # determine read count
                    readcount = linecount // 4 # 0-indexed
                    
                    # loop through each letter in the newline and append phred score to correct list in mean
                    for i, letter in enumerate(newline):
                        phred_score = Bioinfo.convert_phred(letter)
                        mean[i,0] += phred_score

                # increment linecounter
                linecount+=1

    else:

        # open file
        with open(file) as f:
            
            # loop through each line
            while True:
                newline = f.readline().rstrip()
                
                # break if empty (end of file)
                if newline == '':
                    break
                
                # update statement
                if linecount % 500000 == 0:
                    print(linecount)
                
                # grab quality score line of record
                if linecount % 4 == 3:
                    
                    # determine read count
                    readcount = linecount // 4 # 0-indexed
                    
                    # loop through each letter in the newline and add score to running total of mean
                    for i, letter in enumerate(newline):
                        phred_score = Bioinfo.convert_phred(letter)
                        mean[i,0] += phred_score

                # increment linecounter
                linecount+=1
    
    # divide mean array by the number of records to obtain mean
    recordcount: int = int(linecount / 4)
    mean = mean / recordcount

    return mean, linecount

# open file and run array_of_arrays
gzipped: bool = input_file.endswith('.gz')
mean: float, linecount: int = qscore_average_function(
        file=input_file,
        read_length=read_length,
        num_reads=num_reads,
        gzipped=gzipped
    )

# generate and save histogram of mean score by position
plt.bar(range(mean.shape[0]),mean.flatten(),color='purple')
plt.title('mean score by position of read')
plt.xlabel('position')
plt.ylabel('mean score')
plt.savefig(output_file)
