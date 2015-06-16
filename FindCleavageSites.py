# IdentifyOffTargetSiteSequences.py
# Shengdar Tsai (stsai4@mgh.harvard.edu)

# A program to identify Cas9 off-target sites from in vitro GUIDE-seq data

__author__ = 'shengdar'

import argparse
import HTSeq
import itertools
# import pyfaidx

parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
parser.add_argument('--bam', help='Sorted BAM file', required=True)
args = parser.parse_args()

### Tabulate the start positions for the 2nd read in pair across the genome.
### Filter mapped reads with a quality score of 50.

def tabulate_start_positions(BamFileName):
    sorted_bam_file = HTSeq.BAM_Reader(args.bam)
    ga = HTSeq.GenomicArray("auto", stranded=False)
    ga_windows = HTSeq.GenomicArray("auto", stranded=False)
    ga_stranded = HTSeq.GenomicArray("auto", stranded=True)

    for read in itertools.islice( sorted_bam_file, 1000000 ):  # printing first 5 reads
    # for read in sorted_bam_file:
        if read.pe_which == "second" and read.aligned and read.aQual >= 50:
            iv = read.iv
            chr = iv.chrom
            position = iv.start_d
            strand = iv.strand
            ga[HTSeq.GenomicPosition(chr, position, strand)] += 1
            ga_windows[HTSeq.GenomicPosition(chr, position, strand)] = 1
            ga_stranded[HTSeq.GenomicPosition(chr, position, strand)] += 1
    return ga, ga_windows, ga_stranded

### find_windows method
def find_windows(ga, ga_windows, window_size=3):
    # Initialize comparison position
    last = HTSeq.GenomicInterval("0", 0, 0)
    # Iterate through window GenomicArray and consolidate windows that are within 3 bp
    for iv, value in ga_windows.steps():
        if value:
            if iv.chrom != last.chrom or iv.start - last.end > window_size:
                last = iv
            else:
                consolidated_interval = HTSeq.GenomicInterval(iv.chrom, last.start, iv.end)
                ga_windows[consolidated_interval] = 1
                last = consolidated_interval
    # Return consolidated GenomicArray

    for iv, value in ga_windows.steps():
        if value:
            count = sum(list(ga[iv]))
            if count >= 3:
                print iv, sum(list(ga[iv]))

    return ga_windows

def main():
    # Tabulate start positions for read 2
    # Identify positions where there are more than 1 read

    ga, ga_windows, ga_stranded = tabulate_start_positions("/Users/shengdar/Local-projects/invitro_GUIDE-Seq/BAM/SQT01_S1.bam")
    ga_consolidated_windows = find_windows(ga, ga_windows)

if __name__ == "__main__":
    main()


