# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. ```Your answer here```
    3. 
    ```
    Index 1 count of indexes with at least 1 undetermined base:
    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' |  grep N | wc -l
    3976613
    ```
    ```
    Index 2 count of indexes with at least 1 undetermined base:
    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' |  grep N | wc -l
    3328051
    ```


    
## Part 2
1. Define the problem
```
We want to bucket our multiplexed data. We also know that we have some index hopping present in the data, meaning we have to filter out the bad indexes. 
```
2. Describe output
```
The output should include two fastq files for each index-pair: one each for read1 and read2 (24 unique indexes * 2 read files= 48 output fastq files). We should then have two separate fastq files for non-matching index-pairs: one each for read1 and read2 (1 group mismatched * 2 read files = 2 output mismatched index fastq files). We should then output another two fastq files that hold the reads for unknown/low-quality index-pairs: one for each of read1 and read2 (1 group unknown/low-qual * 2 read files = 2 output unknown/low-quality fastq files). Total: 52 output fastq files.

The header for each of the reads in these files must include BOTH index1 and index 2 separated by a hyphen.

We also want to track the number of each of these files. ALSO, for each of the mismatched index, we want to track the number for each possible combination, which is up to 24! number of combinations.
```
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
