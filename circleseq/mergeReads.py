from __future__ import print_function
import argparse
import itertools
import gzip
from utility import reverseComplement, fq

def mergeReads(fastq1_filename, fastq2_filename, out):
    fastq1_file = fq(fastq1_filename)
    fastq2_file = fq(fastq2_filename)

    with gzip.open(out, 'wb') as o:
        for r1, r2 in itertools.izip(fastq1_file, fastq2_file):
            merged_sequence = reverseComplement(r1[1]) + r2[1]
            merged_quality_scores = r1[3][::-1] + r2[3]
            print(r1[0], file=o)
            print(merged_sequence, file=o)
            print(r1[2], file=o)
            print(merged_quality_scores, file=o)

def main():
    parser = argparse.ArgumentParser(description='Merge CIRCLE-seq reads for alignment.')
    parser.add_argument('--read1', help='Read 1 filename', required=True)
    parser.add_argument('--read2', help='Read 2 filename', required=True)
    parser.add_argument('--out', help='Output filename', required=True)

    args = parser.parse_args()

    mergeReads(args.read1, args.read2, args.out)

if __name__ == "__main__":
    main()