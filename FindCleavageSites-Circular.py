# IdentifyOffTargetSiteSequences.py
# Shengdar Tsai (stsai4@mgh.harvard.edu)

# A program to identify Cas9 off-target sites from in vitro GUIDE-seq data


from __future__ import print_function

__author__ = 'shengdar'

import argparse
import HTSeq
import itertools
import os
import pyfaidx
import string
import swalign
import sys

parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
parser.add_argument('--bam', help='Sorted BAM file', required=True)
parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
parser.add_argument('--reads', help='Read threshold', default=2, type=int)
parser.add_argument('--windowsize', help='Read threshold', default=1, type=int)
parser.add_argument('--nofilter', help='Turn off filter for sequence', required=False, action='store_true')

args = parser.parse_args()

### 1. Tabulate the start positions for the 2nd read in pair across the genome.
def tabulate_start_positions(BamFileName):
    sorted_bam_file = HTSeq.BAM_Reader(BamFileName)
    ga = HTSeq.GenomicArray("auto", stranded=False)
    ga_windows = HTSeq.GenomicArray("auto", stranded=False)
    ga_stranded = HTSeq.GenomicArray("auto", stranded=True)
    read_count = 0

    for bundle in HTSeq.pair_SAM_alignments(sorted_bam_file, bundle=True):
        if len(bundle) >= 1:
            for pair in bundle:
                    if pair:
                        first_read, second_read = pair
                        if (first_read is not None) and (second_read is not None):
                             if first_read.aligned and second_read.aligned:
                                first_read_chr = first_read.iv.chrom
                                first_read_position = first_read.iv.start_d
                                first_read_strand = first_read.iv.strand

                                second_read_chr = second_read.iv.chrom
                                second_read_position = second_read.iv.start_d
                                second_read_strand = second_read.iv.strand

                                if (first_read_chr == second_read_chr) and (abs(first_read_position - second_read_position) <= 6):
                                    # print(read_count, first_read_chr, first_read_position, first_read_strand, second_read_chr, second_read_position, second_read_strand)
                                    ga[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1
                                    ga_windows[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] = 1
                                    ga_stranded[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1

                                    ga[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1
                                    ga_windows[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] = 1
                                    ga_stranded[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1

        read_count += 1
        if not read_count % 100000:
            print(read_count/float(1000000), end=" ", file=sys.stderr)

    return ga, ga_windows, ga_stranded

###  2. Find genomic windows (coordinate positions)
def find_windows(ga_windows, window_size):
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

    return ga_windows # Return consolidated GenomicArray

### 3. Find actual sequences of potential off-target sites
def output_alignments(ga, ga_windows, reference_genome):
    target_sequence = args.targetsite
    for iv, value in ga_windows.steps():
        if value:
            count = sum(list(ga[iv]))
            if count >= args.reads:
                window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - 20 , iv.end + 20)
                sequence, mismatches, length, strand,  target_start_relative, target_end_relative = align_sequences(target_sequence, window_sequence)
                if strand == "+":
                    target_start_absolute = target_start_relative + iv.start - 20
                    target_end_absolute = target_end_relative + iv.start - 20
                elif strand == "-":
                    target_start_absolute = iv.end + 20 - target_end_relative
                    target_end_absolute = iv.end + 20 - target_start_relative
                else:
                    target_start_absolute = iv.start
                    target_end_absolute = iv.end
                if sequence or args.nofilter:
                    name = iv.chrom + ':' + str(target_start_absolute) + '-' + str(target_end_absolute)
                    read_count = sum(list(ga[iv]))
                    filename = os.path.basename(args.bam)
                    print(iv.chrom, target_start_absolute, target_end_absolute, name, read_count, strand, iv, iv.chrom,
                          iv.start, iv.end, window_sequence, sequence, mismatches, length, filename, sep="\t")

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

### Get sequences from some reference genome
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][int(start):int(end)]
    elif strand == "-":
        seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
    return str(seq)

### Simple reverse_complement method
def reverse_complement(sequence):
    transtab = string.maketrans("ACGT","TGCA")
    return sequence.translate(transtab)[::-1]

def main():
    # Tabulate start positions for read 2
    # Identify positions where there are more than 1 read
    reference_genome = pyfaidx.Fasta(args.ref)
    print("Reference genome loaded.", file=sys.stderr)
    ga, ga_windows, ga_stranded = tabulate_start_positions(args.bam)
    print("Tabulate start positions.", file=sys.stderr)
    ga_consolidated_windows = find_windows(ga_windows, args.windowsize)
    print("Get consolidated windows.", file=sys.stderr)
    output_alignments(ga, ga_consolidated_windows, reference_genome)
    print("Get alignments.", file=sys.stderr)

if __name__ == "__main__":
    main()


