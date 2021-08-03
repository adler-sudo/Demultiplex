#!/usr/bin/env python

input_file: str = 'barcodes.txt'
output_file: str = 'barcodes.tsv'

with open(input_file,'r') as f, open(output_file,'w') as o:
    line = f.readline()
    while line != '':
        line = line.split()
        for i in range(len(line)):
            if i == 0 or i % 2 == 0:
                o.write(line[i] + "\t")
            else:
                o.write(line[i] + "\n")
        line = f.readline()