# IdentifyOffTargetSiteSequences.py
# Shengdar Tsai (stsai4@mgh.harvard.edu)

# A program to identify Cas9 off-target sites from in vitro GUIDE-seq data

from __future__ import print_function

__author__ = 'shengdar'

import argparse
import HTSeq
import itertools
import pyfaidx
import swalign
import sys

parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
parser.add_argument('--bam', help='Sorted BAM file', required=True)
parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
args = parser.parse_args()

### Tabulate the start positions for the 2nd read in pair across the genome.
### Filter mapped reads with a quality score of 50.
def tabulate_start_positions(BamFileName):
    sorted_bam_file = HTSeq.BAM_Reader(args.bam)
    ga = HTSeq.GenomicArray("auto", stranded=False)
    ga_windows = HTSeq.GenomicArray("auto", stranded=False)
    ga_stranded = HTSeq.GenomicArray("auto", stranded=True)

    for read in itertools.islice( sorted_bam_file, 200000 ):  # printing first N reads
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

### Find genomic windows (coordinate positions)
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

    return ga_windows # Return consolidated GenomicArray

### Find actual sequences of potential off-target sites
def output_alignments(ga, ga_windows, reference_genome):
    for iv, value in ga_windows.steps():
        if value:
            count = sum(list(ga[iv]))
            if count >= 3:
                window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - 25 , iv.end + 25)
                # sequence, mismatches, length, strand,  target_start_relative, target_end_relative = alignSequences(target_sequence, window_sequence)
                print(iv, iv.chrom, iv.start, iv.end, sum(list(ga[iv])), window_sequence, sep="\t" )

### Get sequences from some reference genome
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][int(start):int(end)]
    elif strand == "-":
        seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
    return seq

### Align sequences
def alignSequences(ref_seq, query_seq):
    match = 2
    mismatch = -1
    ref_length = len(ref_seq)
    matches_required = len(ref_seq) - 1 - 7 # allow up to 8 mismatches
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=-100, gap_extension_penalty=-100, prefer_gap_runs=True)  # you can also choose gap penalties, etc...
    forward_alignment = sw.align(ref_seq, query_seq)
    reverse_alignment = sw.align(ref_seq, reverseComplement(query_seq))
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


def main():
    # Tabulate start positions for read 2
    # Identify positions where there are more than 1 read

    reference_genome = pyfaidx.Fasta(args.ref)
    print("Reference genome loaded.", file=sys.stderr)
    ga, ga_windows, ga_stranded = tabulate_start_positions("/Users/shengdar/Local-projects/invitro_GUIDE-Seq/BAM/SQT01_S1.bam")
    print("Tabulate start positions.", file=sys.stderr)
    ga_consolidated_windows = find_windows(ga, ga_windows)
    print("Get consolidated windows.", file=sys.stderr)
    get_alignments(ga, ga_consolidated_windows, reference_genome)
    print("Get alignments.", file=sys.stderr)


if __name__ == "__main__":
    main()


