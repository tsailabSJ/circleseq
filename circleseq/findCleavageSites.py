from __future__ import print_function

__author__ = 'shengdar'

import argparse
import collections
import logging
import HTSeq
import math
import nwalign as nw
import os
import pyfaidx
import regex
import string
import scipy.stats
import numpy
import statsmodels
import sys
import itertools

""" Tabulate merged start positions.
    Identify genomic coordinates for reads mapping across 151/152 bp position.
    Add positions to genomic array.
"""
def tabulate_merged_start_positions(BamFileName, cells, name, targetsite, mapq_threshold, gap_threshold, start_threshold, outfile_base):
    output_filename = '{0}_coordinates.txt'.format(outfile_base)

    sorted_bam_file = HTSeq.BAM_Reader(BamFileName)
    filename_base = os.path.basename(BamFileName)

    ga = HTSeq.GenomicArray("auto", stranded=False)
    ga_windows = HTSeq.GenomicArray("auto", stranded=False)
    ga_stranded = HTSeq.GenomicArray("auto", stranded=True)
    ga_coverage = HTSeq.GenomicArray("auto", stranded=False)

    read_count = 0
    ref_chr = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                '20', '21', '22', 'X', 'Y']


    with open(output_filename, 'w') as o:
        header = ['#Name', 'Targetsite_Sequence', 'Cells', 'BAM', 'Read1_chr', 'Read1_start_position', 'Read1_strand',
                  'Read2_chr', 'Read1_start_position', 'Read2_strand']
        print(*header, sep='\t', file=o)

        for read in sorted_bam_file:
            output = False
            first_read_chr, first_read_position, first_read_strand = None, None, None
            second_read_chr, second_read_position, second_read_strand = None, None, None
            # if not read.flag & 2048 and read.aQual > aqual_threshold:
            if read.aQual > mapq_threshold and read.aligned:

                ga_coverage[read.iv] += 1

                for cigar_operation in read.cigar:
                    # Identify positions that end in position 151 and start at position 151
                    # Note strand polarity is reversed for position 151 (because it is part of the strand that was
                    # reverse complemented initially before alignment
                    if cigar_operation.type == 'M':
                        if (cigar_operation.query_from <= 146 - start_threshold) and (151 - start_threshold <= cigar_operation.query_to):
                            first_read_cigar = cigar_operation
                            first_read_chr = cigar_operation.ref_iv.chrom
                            first_end = min(cigar_operation.query_to, 151)
                            distance = first_end - cigar_operation.query_from
                            first_read_position = cigar_operation.ref_iv.start + distance
                            if cigar_operation.ref_iv.strand == '+':
                                first_read_strand = '-'
                            elif cigar_operation.ref_iv.strand == '-':
                                first_read_strand = '+'
                        if (cigar_operation.query_from <= 151 + start_threshold) and (156 + start_threshold <= cigar_operation.query_to):
                            second_read_cigar = cigar_operation
                            second_read_chr = cigar_operation.ref_iv.chrom
                            second_end = max(151, cigar_operation.query_from)
                            distance = second_end - cigar_operation.query_from
                            if cigar_operation.ref_iv.strand == '+':
                                second_read_position = cigar_operation.ref_iv.start + distance
                                second_read_strand = '+'
                            elif cigar_operation.ref_iv.strand == '-':
                                second_read_position = cigar_operation.ref_iv.start + distance
                                second_read_strand = '-'

                if first_read_chr == second_read_chr and first_read_chr in ref_chr and first_read_position is not None and second_read_position is not None:
                    if abs(first_read_position - second_read_position) <= gap_threshold:
                        output = True
                        ga[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1
                        ga_windows[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] = 1
                        ga_stranded[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1

                        ga[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1
                        ga_windows[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] = 1
                        ga_stranded[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1

                if output == True:
                    print(name, targetsite, cells, filename_base, first_read_chr, first_read_position,
                          first_read_strand, second_read_chr, second_read_position, second_read_strand, sep='\t', file=o)

            read_count += 1
            if not read_count % 100000:
                print(read_count/float(1000000), end=" ", file=sys.stderr)

    return ga, ga_windows, ga_stranded, ga_coverage

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

"""  2. Find genomic windows (coordinate positions)
"""
def find_windows(ga_windows, window_size):
    # Initialize comparison position
    last = HTSeq.GenomicInterval("0", 0, 0)
    # Iterate through window GenomicArray and consolidate windows that are within 3 bp, up to a maximum of 10 bp.
    for iv, value in ga_windows.steps():
        if value:
            if iv.chrom != last.chrom or iv.start - last.end > window_size or iv.end - last.start > 10:
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
    outfile_unmatched = '{0}_identified_unmatched.txt'.format(outfile_base)

    matched_output_table = collections.defaultdict(list)

    #dictionary to store the lines of outfile_matched
    samples_dict = {}   
    #dictionary to add read_count for each pair chromosome:start_position 
    reads_dict = {}
    #dictionary to store window_start. For duplicated off-target sites we will write the minimum of such positions 
    window_min = {}
    #dictionary to store window_end. For duplicated off-target sites we will write the maximum of such positions 
    window_max = {}
    
    with open(outfile_unmatched, 'w') as o2:
        for iv, value in ga_windows.steps():
            if value:
                count = sum(list(ga[iv]))
                if count >= read_threshold:
                    window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - 20 , iv.end + 20)
                    sequence, distance, length, strand,  target_start_relative, target_end_relative = alignSequences(target_sequence, window_sequence, max_errors=6)
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
                        tag = iv.chrom+':'+str(target_start_absolute)
                        if tag not in reads_dict.keys():
                            reads_dict[tag] = read_count
                            window_min[tag] = [iv.start]
                            window_max[tag] = [iv.end]
                            samples_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name, read_count, strand, iv, iv.chrom, iv.start, iv.end,
                                                          window_sequence, sequence, distance, length, filename, target_name, target_cells, full_name, target_sequence]
                        else:
                            reads_dict[tag] = reads_dict[tag] + read_count
                            window_min[tag].append(iv.start)
                            window_max[tag].append(iv.end)
                            samples_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name, reads_dict[tag], strand, iv, iv.chrom, min(window_min[tag]), max(window_max[tag]), 
                                                          window_sequence, sequence, distance, length, filename, target_name, target_cells, full_name, target_sequence]
                    else:
                        print(iv.chrom, target_start_absolute, target_end_absolute, name, read_count, strand, iv, iv.chrom,
                              iv.start, iv.end, window_sequence, sequence, distance, length, filename, target_name,
                              target_cells, full_name,  target_sequence, sep="\t", file=o2)
                        
    return samples_dict



def reverseComplement(sequence):
    transtab = string.maketrans("ACGT","TGCA")
    return sequence.translate(transtab)[::-1]

def regexFromSequence(seq, lookahead=True, indels=1, errors=6):
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
    if errors > 0:
        pattern = pattern + '{{i<={0},d<={0},s<={1},2i+2d+1s<={1}}}'.format(indels, errors)
    return pattern

"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match.
"""
def alignSequences(targetsite_sequence, window_sequence, max_errors=6):
    # Try both strands
    query_regex = regexFromSequence(targetsite_sequence, errors=max_errors)
    forward_alignment = regex.search(query_regex, window_sequence, regex.BESTMATCH)
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
            forward_distance = sum(forward_alignment.fuzzy_counts)
            reverse_distance = sum(reverse_alignment.fuzzy_counts)

            if forward_distance > reverse_distance:
                strand = '-'
                alignment = reverse_alignment
            else:
                strand = '+'
                alignment = forward_alignment

        match_sequence = alignment.group()
        distance = sum(alignment.fuzzy_counts)
        length = len(match_sequence)
        start = alignment.start()
        end = alignment.end()



        if length != len(targetsite_sequence):
            path = os.path.dirname(os.path.abspath(__file__))
            realigned_match_sequence, realigned_target = nw.global_align(match_sequence, targetsite_sequence,
                                                                         gap_open=-10, gap_extend=-100, matrix='{0}/NUC_SIMPLE'.format(path))
            return [realigned_match_sequence, distance, length, strand, start, end]
        else:
            return [match_sequence, distance, length, strand, start, end]

""" Get sequences from some reference genome
"""
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][int(start):int(end)]
    elif strand == "-":
        seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
    return str(seq)


def analyze(ref, bam, targetsite, reads, windowsize, mapq_threshold, gap_threshold, start_threshold, name, cells, out, merged=True):
    output_folder = os.path.dirname(out)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    reference_genome = pyfaidx.Fasta(ref)
    print("Reference genome loaded.", file=sys.stderr)
    print('Reads: {0}, Window: {1}, MAPQ: {2}, Gap: {3}, Start {4}'.format(reads, windowsize, mapq_threshold, gap_threshold, start_threshold), file=sys.stderr)
    if merged:
        print("Tabulate merged start positions.", file=sys.stderr)
        ga, ga_windows, ga_stranded, ga_coverage = tabulate_merged_start_positions(bam, cells, name, targetsite, mapq_threshold, gap_threshold, start_threshold, out)
    else:
        print("Tabulate individual start positions.", file=sys.stderr)
        ga, ga_windows, ga_stranded = tabulate_start_positions(bam, cells, name, targetsite, out)
    ga_consolidated_windows = find_windows(ga_windows, windowsize)
    print("\nGet consolidated windows.", file=sys.stderr)

    samples_dict = output_alignments(ga, ga_consolidated_windows, reference_genome, targetsite, name, cells, bam, reads, out)
    tags_sorted = samples_dict.keys()
    tags_sorted.sort()
    outfile_matched = '{0}_identified_matched.txt'.format(out)
    o1 = open(outfile_matched, 'w')
    for key in tags_sorted:
        print(*samples_dict[key], sep='\t', file=o1)
    o1.close()
    print("Get alignments.", file=sys.stderr)


def compare(ref, bam, control, targetsite, reads, windowsize, mapq_threshold, gap_threshold, start_threshold, name, cells, out, merged=True):

    output_list = list()

    reference_genome = pyfaidx.Fasta(ref)
    combined_ga = HTSeq.GenomicArray("auto", stranded=False) # GenomicArray to store the union of control and nuclease positions
    offtarget_ga_windows = HTSeq.GenomicArray("auto", stranded=False) # GenomicArray to store potential off-target sites

    bg_position = list() # List to store nuclease_position_counts that were observed at least once
    bg_narrow = list() # List to store the sum of nuclease_position_counts in the narrow window
    bg_one_k = list() # List to store the sum of nuclease_position_counts in the one_k window
##    bg_ten_k = list() # List to store the sum of nuclease_position_counts in the ten_k window

    output_folder = os.path.dirname(out)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    output_filename = out + '_count.txt'
    with open(output_filename, 'w') as o:
        if merged:
            print("Tabulate nuclease merged start positions.", file=sys.stderr)
            nuclease_ga, nuclease_ga_windows, nuclease_ga_stranded, nuclease_ga_coverage = tabulate_merged_start_positions(bam,
                cells, name, targetsite, mapq_threshold, gap_threshold, start_threshold, out + '_NUCLEASE')
            print("\nTabulate control merged start positions.", file=sys.stderr)
            control_ga, control_ga_windows, control_ga_stranded, control_ga_coverage = tabulate_merged_start_positions(control,
                cells, name, targetsite, mapq_threshold, gap_threshold, start_threshold, out + '_CONTROL')
            print("\nWriting counts to {0}".format(output_filename), file=sys.stderr)

            # For all positions with detected read mapping positions, put into a combined genomicArray
            for iv, value in nuclease_ga.steps():
                if value:
                    combined_ga[iv] = 1
            for iv, value in control_ga.steps():
                if value:
                    combined_ga[iv] = 1

            for iv, value in combined_ga.steps():
                if value:
                    for position in iv.xrange(step=1):
                        # Define the windows
                        window = HTSeq.GenomicInterval(position.chrom, max(0, position.pos - windowsize), position.pos + windowsize + 1)
                        one_k_window = HTSeq.GenomicInterval(position.chrom, max(0, position.pos - 500), position.pos + windowsize + 500)
                        ten_k_window = HTSeq.GenomicInterval(position.chrom, max(0, position.pos - 5000), position.pos + windowsize + 5000)

                        # Start mapping positions, at the specific base position
                        nuclease_position_counts = nuclease_ga[position]
                        control_position_counts = control_ga[position]
                        # Store control_position_counts for which it was observed at least one read
                        if control_position_counts > 0:
                            bg_position.append(control_position_counts)

                        # In the narrow (parameter-specified) window
                        nuclease_window_counts = sum(nuclease_ga[window])
                        control_window_counts = sum(control_ga[window])
                        # Store control_window_counts greater than zero
                        if control_window_counts > 0:
                            bg_narrow.append(control_window_counts)

                        # In a 1kb window
                        nuclease_one_k_window_counts = sum(nuclease_ga[one_k_window])
                        control_one_k_window_counts = sum(control_ga[one_k_window])
                        # Store control_one_k_window_counts greater than zero
                        if control_one_k_window_counts > 0:
                            bg_one_k.append(control_one_k_window_counts)                       


                        # A list of the outputs, that we will go through again to assign percentiles
                        row = [position.chrom, position.pos, nuclease_position_counts, control_position_counts,
                              nuclease_window_counts, control_window_counts, nuclease_one_k_window_counts, control_one_k_window_counts]
                        output_list.append(row)


            print('#Chromosome', 'zero_based_Position', 'Nuclease_Position_Reads', 'Control_Position_Reads',
                  'Nuclease_Window_Reads', 'Control_Window_Reads',
                  'Nuclease_1k_Window_Reads', 'Control_1k_Window_Reads',
                  'p_Value', 'narrow_p_Value','one_k_p_Value', 'control_p_Value', 'control_narrow_p_Value','control_one_k_p_Value', file=o, sep='\t')


            ecdf_pos = statsmodels.distributions.empirical_distribution.ECDF(bg_position)
            ecdf_nar = statsmodels.distributions.empirical_distribution.ECDF(bg_narrow)
            ecdf_one = statsmodels.distributions.empirical_distribution.ECDF(bg_one_k)
        

            print("\nWritting matched table", file=sys.stderr)

            # Object to store the p-values for every chromosome:position object
            ga_pval = HTSeq.GenomicArray("auto", typecode='O', stranded=False)

            for idx, fields in enumerate(output_list):
                
                position_p_val = 1 - ecdf_pos(fields[2])   
                narrow_p_val = 1 - ecdf_nar(fields[4])     
                one_k_p_val = 1 - ecdf_one(fields[6])      

                control_position_p_val = 1 - ecdf_pos(fields[3])   
                control_narrow_p_val = 1 - ecdf_nar(fields[5])     
                control_one_k_p_val = 1 - ecdf_one(fields[7])    

                if narrow_p_val<0.05 or position_p_val<0.05:
                        
                    read_chr = fields[0]
                    read_position = fields[1]
                    offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
                    
                print(*(fields + [position_p_val, narrow_p_val, one_k_p_val, control_position_p_val, control_narrow_p_val, control_one_k_p_val]), file=o, sep='\t')


                chr_pos = HTSeq.GenomicPosition(fields[0], int(fields[1]), '.')
                ga_pval[chr_pos] = [position_p_val, narrow_p_val, one_k_p_val, control_position_p_val, control_narrow_p_val, control_one_k_p_val]
                                                

            samples_dict = output_alignments(nuclease_ga, offtarget_ga_windows, reference_genome, targetsite, name, cells, bam, reads, out)
            
            tags_sorted = samples_dict.keys()
            tags_sorted.sort()
            outfile_matched = '{0}_identified_matched.txt'.format(out)
            o1 = open(outfile_matched, 'w')
            for key in tags_sorted:

                row = samples_dict[key]       
                
                pos_pval_list = list()
                nar_pval_list = list()
                one_pval_list = list()
                
                control_pos_pval_list = list()
                control_nar_pval_list = list()
                control_one_pval_list = list()                
                 
                iv_pval = HTSeq.GenomicInterval(row[0], int(row[1]), int(row[2]), '.')
                for interval,value in ga_pval[iv_pval].steps():
                    if value is not None:
                        pos_pval_list.append(value[0])
                        nar_pval_list.append(value[1])
                        one_pval_list.append(value[2])
                        control_pos_pval_list.append(value[3])
                        control_nar_pval_list.append(value[4])
                        control_one_pval_list.append(value[5])

                pval_pos = min(pos_pval_list)
                pval_nar = min(nar_pval_list)
                pval_one = min(one_pval_list)
                control_pval_pos = min(control_pos_pval_list)
                control_pval_nar = min(control_nar_pval_list)
                control_pval_one = min(control_one_pval_list)                
                
                    

                print(*(row + [pval_pos, pval_nar, pval_one, control_pval_pos, control_pval_nar, control_pval_one]), sep='\t', file=o1)
            o1.close()
            
            

def main():
    parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
    parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
    parser.add_argument('--bam', help='Sorted BAM file', required=True)
    parser.add_argument('--control', help='Control BAM file', required=False)
    parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
    parser.add_argument('--reads', help='Read threshold', default=4, type=int)
    parser.add_argument('--windowsize', help='Windowsize', default=3, type=int)
    parser.add_argument('--mapq', help='mapq threshold', default=0, type=int)
    parser.add_argument('--gap', help='Gap threshold', default=3, type=int)
    parser.add_argument('--start', help='Start threshold', default=1 , type=int)
    parser.add_argument('--merged', dest='merged', action='store_true', default=False)
    parser.add_argument('--name', help='Targetsite Name', required=False)
    parser.add_argument('--cells', help='Cells', required=False)
    parser.add_argument('--out', help='Output file base', required=True)

    args = parser.parse_args()

    # Run the comparison if the control bam is specified, otherwise run the standard site identification routine.
    if args.control:
        print("Nuclease: {0}\nControl: {1}".format(args.bam, args.control), file=sys.stderr)
        compare(args.ref, args.bam, args.control, args.targetsite, args.reads, args.windowsize, args.mapq, args.gap, args.start, args.name, args.cells, args.out, args.merged)
    else:
        analyze(args.ref, args.bam, args.targetsite, args.reads, args.windowsize, args.mapq, args.gap, args.start, args.name, args.cells, args.out, args.merged)

if __name__ == "__main__":
    main()
