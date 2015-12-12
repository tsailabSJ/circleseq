from __future__ import print_function

import argparse
import itertools
import re
import gzip
import string
import swalign
import collections

"""
FASTQ generator function from umi package
"""
def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]

"""
### Simple reverse_complement method
"""

def reverse_complement(sequence):
    transtab = string.maketrans("ACGT","TGCA")
    return sequence.translate(transtab)[::-1]

"""
### Simple reverse_complement method
"""

### Smith-Waterman alignment of sequences
def align_sequences(ref_seq, query_seq):
    match = 2
    mismatch = -1
    ref_length = len(ref_seq)
    matches_required = len(ref_seq) - 1 - 7 # allow up to 8 mismatches
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=-100, gap_extension_penalty=-100, prefer_gap_runs=True)  # you can also choose gap penalties, etc...
    forward_alignment = sw.align(ref_seq, query_seq)
    reverse_alignment = sw.align(ref_seq, reverse_complement(query_seq))
    if forward_alignment.matches >= matches_required and forward_alignment.matches > reverse_alignment.matches:
        start_pad = forward_alignment.r_pos
        start = forward_alignment.q_pos - start_pad
        end_pad = ref_length - forward_alignment.r_end
        end = forward_alignment.q_end + end_pad
        strand = "+"
        return [forward_alignment.query[start:end], ref_length - forward_alignment.matches - 1, end - start, strand, start, end]
    elif reverse_alignment.matches >= matches_required and reverse_alignment.matches > forward_alignment.matches:
        start_pad = reverse_alignment.r_pos
        start = reverse_alignment.q_pos - start_pad
        end_pad = ref_length - reverse_alignment.r_end
        end = reverse_alignment.q_end + end_pad
        strand = "-"
        return [reverse_alignment.query[start:end], ref_length - reverse_alignment.matches - 1, end - start, strand, start, end]
    else:
        return ["", "", "", "", "", ""]

"""
Main function to find off-target sites in reference-free fashion
"""
def analyze(fastq1_filename, fastq2_filename, targetsite, out, name='', cells=''):

    c = collections.Counter()
    fastq1_file = fq(fastq1_filename)
    fastq2_file = fq(fastq2_filename)
    for r1, r2 in itertools.izip(fastq1_file, fastq2_file):
        r1_sequence = r1[1].rstrip('\n')
        r2_sequence = r2[1].rstrip('\n')
        joined_seq = reverse_complement(r1_sequence) + r2_sequence
        truncated_joined_seq = joined_seq[130:170]
        offtarget, mismatch, length, strand, start, end = align_sequences(targetsite, truncated_joined_seq)
        if offtarget:
            # if strand == '-':
            #     truncated_joined_seq = reverse_complement(truncated_joined_seq)
            # # left_seq = truncated_joined_seq[start - 10:start]
            # # check = truncated_joined_seq[start:end]
            # # right_seq = truncated_joined_seq[end - 10:end]
            c[offtarget] += 1
            print(offtarget, mismatch, length, strand, start, end, c[offtarget])

def join_write_output(fastq1_filename, fastq2_filename, out):
    fastq1_file = fq(fastq1_filename)
    fastq2_file = fq(fastq2_filename)

    with open(out, 'w') as o:
        for r1, r2 in itertools.izip(fastq1_file, fastq2_file):
            header = '>{0}'.format(r1[0])
            r1_sequence = r1[1].rstrip('\n')
            r2_sequence = r2[1].rstrip('\n')
            joined_seq = reverse_complement(r1_sequence) + r2_sequence
            print(header, end='', file=o)
            print(joined_seq, file=o)


def main():
    parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
    parser.add_argument('--fq1', help='FASTQ Read 1', required=True)
    parser.add_argument('--fq2', help='FASTQ Read 2', required=True)
    parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
    parser.add_argument('--name', help='Targetsite Name', required=False)
    parser.add_argument('--cells', help='Cells', required=False)
    parser.add_argument('--out', help='Output file base', required=True)
    args = parser.parse_args()

    analyze(args.fq1, args.fq2, args.targetsite, args.out, args.name, args.cells)
    # join_write_output(args.fq1, args.fq2, args.out)

if __name__ == "__main__":
    main()
