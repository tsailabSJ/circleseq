from __future__ import print_function

__author__ = 'shengdar'

import argparse
import collections
import HTSeq
import os
import pyfaidx
import regex
import string
import sys


""" Tabulate the start positions for the 2nd read in pair across the genome.
    Only consider alignments with matching positions from the beginning of the read.
    For read pairs with multiple alignments, pick the one with matching positions at the beginning.
"""
def tabulate_start_positions(BamFileName, cells, name, targetsite, outfile_base):

    output_filename = '{0}_coordinates.txt'.format(outfile_base)

    sorted_bam_file = HTSeq.BAM_Reader(BamFileName)
    filename_base = os.path.basename(BamFileName)
    ga = HTSeq.GenomicArray("auto", stranded=False)
    ga_windows = HTSeq.GenomicArray("auto", stranded=False)
    ga_stranded = HTSeq.GenomicArray("auto", stranded=True)
    read_count = 0
    current_pair_position = []
    last_pair_position = []
    ref_chr = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                '20', '21', '22', 'X', 'Y']
    aqual_threshold = 0

    with open(output_filename, 'w') as o:
        header = ['#Name', 'Targetsite_Sequence', 'Cells', 'BAM', 'Read1_chr', 'Read1_start_position', 'Read1_strand',
                  'Read2_chr', 'Read1_start_position', 'Read2_strand']
        print(*header, sep='\t', file=o)
        for bundle in HTSeq.pair_SAM_alignments(sorted_bam_file, bundle=True):
            output = False
            first_read_chr, first_read_position, first_read_strand = None, None, None
            second_read_chr, second_read_position, second_read_strand = None, None, None

            if len(bundle) == 1: # single alignment
                first_read, second_read = bundle[0]
                if first_read.aligned:
                    if first_read.aQual >= aqual_threshold and not first_read.flag & 1024 and \
                    (first_read.iv.strand == '+' and first_read.cigar[0].type == 'M') or \
                    (first_read.iv.strand == '-' and first_read.cigar[-1].type == 'M'):
                        first_read_chr = first_read.iv.chrom
                        first_read_position = first_read.iv.start_d
                        first_read_strand = first_read.iv.strand
                if second_read.aligned:
                    if second_read.aQual >= aqual_threshold and not first_read.flag & 1024 and \
                    (second_read.iv.strand == '+' and second_read.cigar[0].type == 'M') or \
                    (second_read.iv.strand == '-' and second_read.cigar[-1].type == 'M'):
                        second_read_chr = second_read.iv.chrom
                        second_read_position = second_read.iv.start_d
                        second_read_strand = second_read.iv.strand
            elif len(bundle) > 1: # multiple alignments
                first_read_list, second_read_list = zip(*bundle)
                filtered_first_read_list = []
                filtered_second_read_list = []
                for read in first_read_list:
                    if read:
                        if read.aligned:
                            if read.iv.strand == '+' and read.cigar[0].type == 'M':
                                    filtered_first_read_list.append(read)
                            elif read.iv.strand == '-' and read.cigar[-1].type == 'M':
                                    filtered_first_read_list.append(read)
                for read in second_read_list:
                    if read:
                        if read.aligned:
                            if read.iv.strand == '+' and read.cigar[0].type == 'M':
                                    filtered_second_read_list.append(read)
                            elif read.iv.strand == '-' and read.cigar[-1].type == 'M':
                                    filtered_second_read_list.append(read)
                if len(filtered_first_read_list) == 1:
                    first_read = filtered_first_read_list[0]
                    if first_read.aQual >= aqual_threshold and not first_read.flag & 1024:
                        first_read_chr = first_read.iv.chrom
                        first_read_position = first_read.iv.start_d
                        first_read_strand = first_read.iv.strand
                if len(filtered_second_read_list) == 1:
                    second_read = filtered_second_read_list[0]
                    if second_read.aQual >= aqual_threshold and not first_read.flag & 1024:
                        second_read_chr = second_read.iv.chrom
                        second_read_position = second_read.iv.start_d
                        second_read_strand = second_read.iv.strand

            # Only count pairs where they are on the same chromosome, originate from within 6 bp start positions,


            # current_pair_position = [first_read_chr, first_read_position, first_read_strand, second_read_chr, second_read_position, second_read_strand]
            # if first_read_chr == second_read_chr and first_read_chr in ref_chr and current_pair_position != last_pair_position and \
            # ((first_read.iv.strand == '+' and second_read.iv.strand == '-' and abs(first_read_position - second_read_position - 1) <= 20)
            # or (second_read.iv.strand == '+' and first_read.iv.strand == '-' and abs(second_read_position - first_read_position - 1) <= 20)):
            # if current_pair_position != last_pair_position:

            if first_read_chr == second_read_chr and first_read_chr in ref_chr and \
            ((first_read.iv.strand == '+' and second_read.iv.strand == '-' and abs(first_read_position - second_read_position - 1) <= 20)
            or (second_read.iv.strand == '+' and first_read.iv.strand == '-' and abs(second_read_position - first_read_position - 1) <= 20)):

                #if first_read_chr in ref_chr and first_read_position and first_read_strand:
                ga[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1
                ga_windows[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] = 1
                ga_stranded[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1
                #    output = True

                #if second_read_chr in ref_chr and second_read_position and second_read_strand:
                ga[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1
                ga_windows[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] = 1
                ga_stranded[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1
                output = True

            # Output read positions for plotting. Add gap.

            if output == True:
                print(name, targetsite, cells, filename_base, first_read_chr, first_read_position,
                      first_read_strand, second_read_chr, second_read_position, second_read_strand, sep='\t', file=o)

            last_pair_position = [ first_read_chr, first_read_position, first_read_strand, second_read_chr, second_read_position, second_read_strand]

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

""" 3. Find actual sequences of potential off-target sites
"""
def output_alignments(ga, ga_windows, reference_genome, target_sequence, target_name, target_cells, bam_filename, read_threshold, outfile_base):

    # outfile_dirname, outfile_basename = os.path.split(outfile_base)
    outfile_matched = '{0}_identified_matched.txt'.format(outfile_base)
    outfile_unmatched = '{0}_identified_unmatched.txt'.format(outfile_base)

    matched_output_table = collections.defaultdict(list)

    with open(outfile_matched, 'w') as o1, open(outfile_unmatched, 'w') as o2:
        for iv, value in ga_windows.steps():
            if value:
                count = sum(list(ga[iv]))
                if count >= read_threshold:
                    window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - 20 , iv.end + 20)
                    sequence, mismatches, length, strand,  target_start_relative, target_end_relative = alignSequences(target_sequence, window_sequence, 7)
                    if strand == "+":
                        target_start_absolute = target_start_relative + iv.start - 20
                        target_end_absolute = target_end_relative + iv.start - 20
                    elif strand == "-":
                        target_start_absolute = iv.end + 20 - target_end_relative
                        target_end_absolute = iv.end + 20 - target_start_relative
                    else:
                        target_start_absolute = iv.start
                        target_end_absolute = iv.end
                    name = iv.chrom + ':' + str(target_start_absolute) + '-' + str(target_end_absolute)
                    read_count = int(sum(list(ga[iv])))
                    filename = os.path.basename(bam_filename)
                    full_name = target_name + '_' + target_cells + '_' + name + '_' + str(read_count)
                    if sequence:
                        print(iv.chrom, target_start_absolute, target_end_absolute, name, read_count, strand, iv, iv.chrom,
                              iv.start, iv.end, window_sequence, sequence, mismatches, length, filename, target_name,
                              target_cells, full_name, target_sequence, sep="\t", file=o1)
                    else:
                        print(iv.chrom, target_start_absolute, target_end_absolute, name, read_count, strand, iv, iv.chrom,
                              iv.start, iv.end, window_sequence, sequence, mismatches, length, filename, target_name,
                              target_cells, full_name,  target_sequence, sep="\t", file=o2)

def reverseComplement(sequence):
    transtab = string.maketrans("ACGT","TGCA")
    return sequence.translate(transtab)[::-1]

def regexFromSequence(seq, lookahead=True, indels=1, mismatches=2):
    """
    Given a sequence with ambiguous base characters, returns a regex that matches for
    the explicit (unambiguous) base characters
    """
    IUPAC_notation_regex = {'N': '[ATCGN]',
                            'Y': '[CTY]',
                            'R': '[AGR]',
                            'W': '[ATW]',
                            'S': '[CGS]',
                            'A': 'A',
                            'T': 'T',
                            'C': 'C',
                            'G': 'G'}

    pattern = ''

    for c in seq:
        pattern += IUPAC_notation_regex[c]

    if lookahead:
        pattern = '(?:' + pattern + ')'
    if mismatches > 0:
        pattern = pattern + '{{s<={}}}'.format(mismatches)
        # pattern = pattern + '{{i<={0},d<={1},1i+1d+s<={2}}}'.format(indels, indels, mismatches)
    return pattern

"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match.
"""
def alignSequences(targetsite_sequence, window_sequence, max_mismatches = 6):
    # Try both strands
    query_regex = regexFromSequence(targetsite_sequence, mismatches=max_mismatches)
    forward_alignment = regex.search(query_regex, window_sequence, regex.BESTMATCH)

    # reverse_regex = regexFromSequence(reverseComplement(targetsite_sequence), mismatches=max_mismatches)
    reverse_alignment = regex.search(query_regex, reverseComplement(window_sequence), regex.BESTMATCH)

    if forward_alignment is None and reverse_alignment is None:
        return ['', '', '', '', '', '']
    else:
        if forward_alignment is None and reverse_alignment is not None:
            strand = '-'
            alignment = reverse_alignment
        elif reverse_alignment is None and forward_alignment is not None:
            strand = '+'
            alignment = forward_alignment
        elif forward_alignment is not None and reverse_alignment is not None:
            forward_mismatches = forward_alignment.fuzzy_counts[0]
            reverse_mismatches = reverse_alignment.fuzzy_counts[0]

            if forward_mismatches > reverse_mismatches:
                strand = '-'
                alignment = reverse_alignment
            else:
                strand = '+'
                alignment = forward_alignment

        match_sequence = alignment.group()
        mismatches = alignment.fuzzy_counts[0]
        length = len(match_sequence)
        start = alignment.start()
        end = alignment.end()

        return [match_sequence, mismatches, length, strand, start, end]

### Get sequences from some reference genome
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][int(start):int(end)]
    elif strand == "-":
        seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
    return str(seq)

def analyze(ref, bam, targetsite, reads, windowsize, name, cells, out, merged=False):
    output_folder = os.path.dirname(out)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    reference_genome = pyfaidx.Fasta(ref)
    print("Reference genome loaded.", file=sys.stderr)
    if merged:
        pass
    else:
        ga, ga_windows, ga_stranded = tabulate_start_positions(bam, cells, name, targetsite, out)
        print("Tabulate start positions.", file=sys.stderr)
        ga_consolidated_windows = find_windows(ga_windows, windowsize)
        print("Get consolidated windows.", file=sys.stderr)
        output_alignments(ga, ga_consolidated_windows, reference_genome, targetsite, name, cells, bam, reads, out)
        print("Get alignments.", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
    parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
    parser.add_argument('--bam', help='Sorted BAM file', required=True)
    parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
    parser.add_argument('--reads', help='Read threshold', default=4, type=int)
    parser.add_argument('--windowsize', help='Windowsize', default=3, type=int)
    parser.add_argument('--merged', dest='merged', action='store_true', default=False)
    # parser.add_argument('--nofilter', help='Turn off filter for sequence', required=False, action='store_true')
    parser.add_argument('--name', help='Targetsite Name', required=False)
    parser.add_argument('--cells', help='Cells', required=False)
    parser.add_argument('--out', help='Output file base', required=True)

    args = parser.parse_args()

    analyze(args.ref, args.bam, args.targetsite, args.reads, args.windowsize, args.name, args.cells, args.out, args.merged)

if __name__ == "__main__":
    main()
