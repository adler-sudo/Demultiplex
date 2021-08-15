# Demultiplexing and Index Swapping - Assignment the Third

Write your code to demultiplex the samples. Be sure to:

- Incorporate feedback from peer code reviews
- Utilize appropriate functions
- Sufficiently comment your code/use docstrings
- Use unit tests on functions/entire algorithm to ensure it works properly
- Create a useful report for the end user of your code
- Use argparse to "generalize" your code
- Follow the specifications laid out in [Assignment the First](../Assignment-the-first) for the code
    - Unclear? Ask Leslie

Final work will be submitted on [GitHub](.). Make sure your folder is well organized and final output is clearly labeled/summarized (a markdown file would be much appreciated!!). Use your code to demultiplex the samples and report:
- Percentage of reads from each sample
- Overall amount of index swapping
- Any figures/any other relevant data your code output

# Description of items
Bioinfo.py - Bioinfo module
barcode_parser.py - generates tsv from barcodes.txt
barcodes.tsv - tsv of barcodes for current multiplexed run
barcodes.txt - txt of barcodes for current multiplexed run
demux.out - output from talapas
demux.py - demultiplex algorithm
demux.srun - the instructions sent to talapas
summary_stats.md - summary output of the current multiplexed run
