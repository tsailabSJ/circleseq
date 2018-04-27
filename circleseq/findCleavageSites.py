from __future__ import print_function

import argparse
import HTSeq
import os
import pyfaidx
import regex
from statsmodels.distributions.empirical_distribution import ECDF
import sys


""" Tabulate merged start positions.
    Identify genomic coordinates for reads mapping across 151/152 bp position.
    Add positions to genomic array.
"""
def tabulate_merged_start_positions(BamFileName, cells, name, targetsite, mapq_threshold, gap_threshold, start_threshold, outfile_base, pattern, all_chromosomes):
    output_filename = '{0}_coordinates.txt'.format(outfile_base)

    sorted_bam_file = HTSeq.BAM_Reader(BamFileName)
    filename_base = os.path.basename(BamFileName)

    ga = HTSeq.GenomicArray("auto", stranded=False)
    ga_windows = HTSeq.GenomicArray("auto", stranded=False)
    ga_stranded = HTSeq.GenomicArray("auto", stranded=True)
    ga_coverage = HTSeq.GenomicArray("auto", stranded=False)

    read_count = 0

    with open(output_filename, 'w') as o:
        header = ['#Name', 'Targetsite_Sequence', 'Cells', 'BAM', 'Read1_chr', 'Read1_start_position', 'Read1_strand',
                  'Read2_chr', 'Read1_start_position', 'Read2_strand']
        print(*header, sep='\t', file=o)

        for read in sorted_bam_file:
            output = False
            first_read_chr, first_read_position, first_read_strand = None, None, None
            second_read_chr, second_read_position, second_read_strand = None, None, None
            if read.aQual > mapq_threshold and read.aligned:

                ga_coverage[read.iv] += 1

                for cigar_operation in read.cigar:
                    # Identify positions that end in position 151 and start at position 151
                    # Note strand polarity is reversed for position 151 (because it is part of the strand that was
                    # reverse complemented initially before alignment
                    if cigar_operation.type == 'M':
                        if ((cigar_operation.query_from <= 146 - start_threshold) and
                                (151 - start_threshold <= cigar_operation.query_to)):
                            first_read_cigar = cigar_operation
                            first_read_chr = cigar_operation.ref_iv.chrom
                            first_end = min(cigar_operation.query_to, 151)
                            distance = first_end - cigar_operation.query_from
                            first_read_position = cigar_operation.ref_iv.start + distance - 1
                            first_read_strand = '-'
                        if ((cigar_operation.query_from <= 151 + start_threshold) and
                                (156 + start_threshold <= cigar_operation.query_to)):
                            second_read_cigar = cigar_operation
                            second_read_chr = cigar_operation.ref_iv.chrom
                            second_end = max(151, cigar_operation.query_from)
                            distance = second_end - cigar_operation.query_from
                            second_read_position = cigar_operation.ref_iv.start + distance
                            second_read_strand = '+'

                if first_read_chr:
                    if first_read_chr == second_read_chr and (pattern.match(str(first_read_chr)) or all_chromosomes) and \
                                    first_read_position is not None and second_read_position is not None:
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

    return ga, ga_windows, ga_stranded, ga_coverage, read_count


""" Tabulate the start positions for the 2nd read in pair across the genome.
    Only consider alignments with matching positions from the beginning of the read.
    For read pairs with multiple alignments, pick the one with matching positions at the beginning.
"""
def tabulate_start_positions(BamFileName, cells, name, targetsite, mapq_threshold, gap_threshold, outfile_base, pattern, all_chromosomes):

    output_filename = '{0}_coordinates.txt'.format(outfile_base)

    sorted_bam_file = HTSeq.BAM_Reader(BamFileName)
    filename_base = os.path.basename(BamFileName)
    ga = HTSeq.GenomicArray("auto", stranded=False)
    ga_windows = HTSeq.GenomicArray("auto", stranded=False)
    ga_stranded = HTSeq.GenomicArray("auto", stranded=True)
    ga_coverage = HTSeq.GenomicArray("auto", stranded=False)
    read_count = 0

    with open(output_filename, 'w') as o:
        header = ['#Name', 'Targetsite_Sequence', 'Cells', 'BAM', 'Read1_chr', 'Read1_start_position', 'Read1_strand',
                  'Read2_chr', 'Read1_start_position', 'Read2_strand']
        print(*header, sep='\t', file=o)
        for bundle in HTSeq.pair_SAM_alignments(sorted_bam_file, bundle=True):
            output = False
            first_read_chr, first_read_position, first_read_strand = None, None, None
            second_read_chr, second_read_position, second_read_strand = None, None, None

            if len(bundle) == 1:  # single alignment
                first_read, second_read = bundle[0]
                if first_read.aligned:
                    if first_read.aQual >= mapq_threshold and not first_read.flag & 1024 and \
                            ((first_read.iv.strand == '+' and first_read.cigar[0].type == 'M') or (first_read.iv.strand == '-' and first_read.cigar[-1].type == 'M')):
                        first_read_chr = first_read.iv.chrom
                        first_read_position = first_read.iv.start_d
                        first_read_strand = first_read.iv.strand
                if second_read.aligned:
                    if second_read.aQual >= mapq_threshold and not first_read.flag & 1024 and \
                            ((second_read.iv.strand == '+' and second_read.cigar[0].type == 'M') or (second_read.iv.strand == '-' and second_read.cigar[-1].type == 'M')):
                        second_read_chr = second_read.iv.chrom
                        second_read_position = second_read.iv.start_d
                        second_read_strand = second_read.iv.strand
            elif len(bundle) > 1:  # multiple alignments
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
                    if first_read.aQual >= mapq_threshold and not first_read.flag & 1024:
                        first_read_chr = first_read.iv.chrom
                        first_read_position = first_read.iv.start_d
                        first_read_strand = first_read.iv.strand
                if len(filtered_second_read_list) == 1:
                    second_read = filtered_second_read_list[0]
                    if second_read.aQual >= mapq_threshold and not first_read.flag & 1024:
                        second_read_chr = second_read.iv.chrom
                        second_read_position = second_read.iv.start_d
                        second_read_strand = second_read.iv.strand

            # We check whether or not the read was aligned by asking for 'first_read_chr'
            if first_read_chr:
                if (first_read_chr == second_read_chr) and (pattern.match(str(first_read_chr)) or all_chromosomes) and \
                        ((first_read.iv.strand == '+' and second_read.iv.strand == '-' and abs(first_read_position - second_read_position) <= gap_threshold) or
                             (second_read.iv.strand == '+' and first_read.iv.strand == '-' and abs(second_read_position - first_read_position) <= gap_threshold)):

                    # if first_read_chr in ref_chr and first_read_position and first_read_strand:
                    ga[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1
                    ga_windows[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] = 1
                    ga_stranded[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1
                    #    output = True

                    # if second_read_chr in ref_chr and second_read_position and second_read_strand:
                    ga[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1
                    ga_windows[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] = 1
                    ga_stranded[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1
                    output = True

            # Output read positions for plotting. Add gap.
            if output == True:
                print(name, targetsite, cells, filename_base, first_read_chr, first_read_position,
                      first_read_strand, second_read_chr, second_read_position, second_read_strand, sep='\t', file=o)

            last_pair_position = [first_read_chr, first_read_position, first_read_strand, second_read_chr, second_read_position, second_read_strand]

            read_count += 1
            if not read_count % 100000:
                print(read_count/float(1000000), end=" ", file=sys.stderr)

    return ga, ga_windows, ga_stranded, ga_coverage, read_count

""" Find genomic windows (coordinate positions)
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

    return ga_windows  # Return consolidated GenomicArray

""" Find actual sequences of potential off-target sites
"""
def output_alignments(narrow_ga, ga_windows, reference_genome, target_sequence, target_name, target_cells,
                      bam_filename, mismatch_threshold, ga_pval, search_radius, out):

    # dictionary to store the matched reads
    matched_dict = {}   
    # dictionary to add read_count for each pair chromosome:start_position among matched reads
    reads_dict = {}
    # dictionary to store window_start. For duplicated matched off-target.
    window_min = {}
    # dictionary to store window_end. For duplicated matched off-target.
    window_max = {}

    # dictionary to store the unmatched reads
    unmatched_dict = {}
    
    for iv, value in ga_windows.steps():
        if value:
            window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - search_radius, iv.end + search_radius)

            offtarget_sequence_no_bulge, mismatches, offtarget_sequence_length, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
            realigned_target, \
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end = \
                alignSequences(target_sequence, window_sequence, max_score=mismatch_threshold)

            # get genomic coordinates of sequences
            mm_start, mm_end, b_start, b_end = '', '', '', ''
            if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '+':
                mm_start = iv.start - search_radius + int(start_no_bulge)
                mm_end = iv.start - search_radius + int(end_no_bulge)
            if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '-':
                mm_start = iv.end + search_radius - int(end_no_bulge)
                mm_end = iv.end + search_radius - int(start_no_bulge)

            if bulged_offtarget_sequence and chosen_alignment_strand_b == '+':
                b_start = iv.start - search_radius + int(bulged_start)
                b_end = iv.start - search_radius + int(bulged_end)
            if bulged_offtarget_sequence and chosen_alignment_strand_b == '-':
                b_start = iv.end + search_radius - int(bulged_end)
                b_end = iv.end + search_radius - int(bulged_start)

            #  define overall start, end and strand. For bed annotation purposes
            if offtarget_sequence_no_bulge:
                target_start_absolute = mm_start
                target_end_absolute = mm_end
                target_strand_absolute = chosen_alignment_strand_m
            elif not offtarget_sequence_no_bulge and bulged_offtarget_sequence:
                target_start_absolute = b_start
                target_end_absolute = b_end
                target_strand_absolute = chosen_alignment_strand_b
            else:
                target_start_absolute = iv.start
                target_end_absolute = iv.end
                target_strand_absolute = '*'

            name = iv.chrom + ':' + str(target_start_absolute) + '-' + str(target_end_absolute)
            read_count = int(max(set(narrow_ga[iv])))
            filename = os.path.basename(bam_filename)
            full_name = str(target_name) + '_' + str(target_cells) + '_' + str(name) + '_' + str(read_count)

            if offtarget_sequence_no_bulge or bulged_offtarget_sequence:
                tag = iv.chrom + ':' + str(target_start_absolute)
                if tag not in reads_dict.keys():
                    reads_dict[tag] = read_count
                    window_min[tag] = [iv.start]
                    window_max[tag] = [iv.end]
                    matched_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name, read_count, target_strand_absolute,
                                         iv.start, iv.end, iv, window_sequence,
                                         offtarget_sequence_no_bulge, mismatches,
                                         chosen_alignment_strand_m, mm_start, mm_end,
                                         bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
                                         chosen_alignment_strand_b, b_start, b_end,
                                         filename, target_cells, target_name, full_name, target_sequence, realigned_target]
                else:
                    current_read_count = reads_dict[tag]
                    reads_dict[tag] = max(current_read_count, read_count) 
                    window_min[tag].append(iv.start)
                    window_max[tag].append(iv.end)
                    matched_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name, read_count, target_strand_absolute,
                                         min(window_min[tag]), max(window_max[tag]), iv, window_sequence,
                                         offtarget_sequence_no_bulge, mismatches,
                                         chosen_alignment_strand_m, mm_start, mm_end,
                                         bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
                                         chosen_alignment_strand_b, b_start, b_end,
                                         filename, target_cells, target_name, full_name, target_sequence, realigned_target]
            else:
                untag = iv.chrom + ':' + str(iv.start)
                unmatched_dict[untag] = [iv.chrom, target_start_absolute, target_end_absolute, name, read_count, target_strand_absolute,
                                         iv.start, iv.end, iv, window_sequence,
                                         offtarget_sequence_no_bulge, mismatches,
                                         chosen_alignment_strand_m, mm_start, mm_end,
                                         bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
                                         chosen_alignment_strand_b, b_start, b_end,
                                         filename, target_cells, target_name, full_name, target_sequence, 'none']

    # Write matched table
    print("Writing matched table", file=sys.stderr)
    tags_sorted = matched_dict.keys()
    tags_sorted.sort()
    outfile_matched = '{0}_identified_matched.txt'.format(out)

    o1 = open(outfile_matched, 'w')
    print('Chromosome', 'Start', 'End', 'Name', 'ReadCount', 'Strand',  # 0:5
          'MappingPositionStart', 'MappingPositionEnd', 'WindowName', 'WindowSequence',  # 6:9
          'Site_SubstitutionsOnly.Sequence', 'Site_SubstitutionsOnly.NumSubstitutions',  # 10:11
          'Site_SubstitutionsOnly.Strand', 'Site_SubstitutionsOnly.Start', 'Site_SubstitutionsOnly.End',  # 12:14
          'Site_GapsAllowed.Sequence', 'Site_GapsAllowed.Length', 'Site_GapsAllowed.Score',  # 15:17
          'Site_GapsAllowed.Substitutions', 'Site_GapsAllowed.Insertions', 'Site_GapsAllowed.Deletions',  # 18:20
          'Site_GapsAllowed.Strand', 'Site_GapsAllowed.Start', 'Site_GapsAllowed.End',  #21:23
          'FileName', 'Cell', 'Targetsite', 'FullName', 'TargetSequence', 'RealignedTargetSequence',  # 24:29
          'Position.Pvalue', 'Narrow.Pvalue', 'Position.Control.Pvalue', 'Narrow.Control.Pvalue',  # 30:33
          sep='\t', file=o1)
    o1.close()

    with open(outfile_matched, 'a') as o1:
        for key in tags_sorted:
            row = matched_dict[key]

            pos_pval_list, nar_pval_list = list(), list()
            control_pos_pval_list, control_nar_pval_list = list(), list()

            iv_pval = HTSeq.GenomicInterval(row[0], int(row[1]), int(row[2]), '.')
            for interval, value in ga_pval[iv_pval].steps():
                if value is not None:
                    pos_pval_list.append(value[0])
                    nar_pval_list.append(value[1])
                    control_pos_pval_list.append(value[2])
                    control_nar_pval_list.append(value[3])

            pval_pos = min(pos_pval_list)
            pval_nar = min(nar_pval_list)
            control_pval_pos = min(control_pos_pval_list)
            control_pval_nar = min(control_nar_pval_list)

            print(*(row + [pval_pos, pval_nar, control_pval_pos, control_pval_nar]), sep='\t', file=o1)

    # Write unmatched table
    print("Writing unmatched table", file=sys.stderr)
    untags_sorted = unmatched_dict.keys()
    untags_sorted.sort()
    outfile_unmatched = '{0}_identified_unmatched.txt'.format(out)
    with open(outfile_unmatched, 'w') as o2:
        for unkey in untags_sorted:
            unrow = unmatched_dict[unkey]

            un_pos_pval_list, un_nar_pval_list = list(), list()
            un_control_pos_pval_list, un_control_nar_pval_list = list(), list()

            iv_pval = HTSeq.GenomicInterval(unrow[0], int(unrow[1]), int(unrow[2]), '.')
            for interval, value in ga_pval[iv_pval].steps():
                if value is not None:
                    un_pos_pval_list.append(value[0])
                    un_nar_pval_list.append(value[1])
                    un_control_pos_pval_list.append(value[2])
                    un_control_nar_pval_list.append(value[3])

            un_pval_pos = min(un_pos_pval_list)
            un_pval_nar = min(un_nar_pval_list)
            un_control_pval_pos = min(un_control_pos_pval_list)
            un_control_pval_nar = min(un_control_nar_pval_list)

            print(*(unrow + [un_pval_pos, un_pval_nar, un_control_pval_pos, un_control_pval_nar]), sep='\t', file=o2)


""" Reverse complement DNA sequence
"""
def reverseComplement(seq):
    compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', '.': '.', '-': '-', '_': '_'})
    out_list = [compl[bp] for bp in seq]
    return ''.join(out_list[::-1])


def regexFromSequence(seq, lookahead=True, indels=1, errors=7):
    seq = seq.upper()
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
        pattern = '(?b:' + pattern + ')'

    pattern_standard = pattern + '{{s<={0}}}'.format(errors)
    pattern_gap = pattern + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)
    return pattern_standard, pattern_gap

"""
Allow for '-' in our search, but do not allow insertions or deletions. 
"""
def extendedPattern(seq, errors):
    IUPAC_notation_regex_extended = {'N': '[ATCGN]','-': '[ATCGN]','Y': '[CTY]','R': '[AGR]','W': '[ATW]','S': '[CGS]','A': 'A','T': 'T','C': 'C','G': 'G'}
    realign_pattern = ''
    for c in seq:
        realign_pattern += IUPAC_notation_regex_extended[c]
    return '(?b:' + realign_pattern + ')' + '{{s<={0}}}'.format(errors)

"""
Recreate A!!! sequence in the window_sequence that matches the conditions given for the fuzzy regex. 
Currently only working for off-targets with at most one bulge !!! 
"""
def realignedSequences(targetsite_sequence, chosen_alignment, errors):
    match_sequence = chosen_alignment.group()
    substitutions, insertions, deletions = chosen_alignment.fuzzy_counts

    # get the .fuzzy_counts associated to the matching sequence after adjusting for indels, where 0 <= INS, DEL <= 1
    realigned_fuzzy = (substitutions, max(0, insertions - 1), max(0, deletions - 1))

    #  recreate the target sequence, with a '-' in the case of an DNA-bulge
    if insertions:
        targetsite_realignments = [targetsite_sequence[:i] + '-' + targetsite_sequence[i:] for i in range(1, len(targetsite_sequence))]
    else:
        targetsite_realignments = [targetsite_sequence]

    realigned_target_sequence, realigned_offtarget_sequence = None, ''  # in case the matching sequence is not founded

    for seq in targetsite_realignments:
        # recreate off-target sequence (with a '-' in the case of an RNA-bulge) and pattern matching the realigned target sequence
        if deletions:
            match_realignments = [match_sequence[:i + 1] + '-' + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
            match_pattern = [match_sequence[:i + 1] + seq[i + 1] + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
        else:
            match_realignments = match_pattern = [match_sequence]

        x = extendedPattern(seq, errors)
        for y_pattern, y_alignment in zip(match_pattern, match_realignments):
            m = regex.search(x, y_pattern, regex.BESTMATCH)
            if m and m.fuzzy_counts == realigned_fuzzy:
                realigned_target_sequence, realigned_offtarget_sequence = seq, y_alignment
    return realigned_target_sequence, realigned_offtarget_sequence


"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match(es).
"""
def alignSequences(targetsite_sequence, window_sequence, max_score=7):

    window_sequence = window_sequence.upper()
    query_regex_standard, query_regex_gap = regexFromSequence(targetsite_sequence, errors=max_score)

    # Try both strands
    alignments_mm, alignments_bulge = list(), list()
    alignments_mm.append(('+', 'standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
    alignments_mm.append(('-', 'standard', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
    alignments_bulge.append(('+', 'gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
    alignments_bulge.append(('-', 'gapped', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))

    lowest_distance_score, lowest_mismatch = 100, max_score + 1
    chosen_alignment_b, chosen_alignment_m, chosen_alignment_strand_b, chosen_alignment_strand_m = None, None, '', ''

    # Use regex to find the best match allowing only for mismatches
    for aln_m in alignments_mm:
        strand_m, alignment_type_m, match_m = aln_m
        if match_m != None:
            mismatches, insertions, deletions = match_m.fuzzy_counts
            if mismatches < lowest_mismatch:
                chosen_alignment_m = match_m
                chosen_alignment_strand_m = strand_m
                lowest_mismatch = mismatches

    # Use regex to find the best match allowing for gaps, so that its edit distance is strictly lower than the
    # total number of mismatches of the sequence founded (if any) allowing only for mismatches.
    for aln_b in alignments_bulge:
        strand_b, alignment_type_b, match_b = aln_b
        if match_b != None:
            substitutions, insertions, deletions = match_b.fuzzy_counts
            if insertions or deletions:
                distance_score = substitutions + (insertions + deletions) * 3
                edistance = substitutions + insertions + deletions
                if distance_score < lowest_distance_score and edistance < lowest_mismatch:
                    chosen_alignment_b = match_b
                    chosen_alignment_strand_b = strand_b
                    lowest_distance_score = distance_score

    if chosen_alignment_m:
        offtarget_sequence_no_bulge = chosen_alignment_m.group()
        mismatches = chosen_alignment_m.fuzzy_counts[0]
        start_no_bulge = chosen_alignment_m.start()
        end_no_bulge = chosen_alignment_m.end()
    else:
        offtarget_sequence_no_bulge, mismatches, start_no_bulge, end_no_bulge, chosen_alignment_strand_m = '', '', '', '', ''

    bulged_offtarget_sequence, score, length, substitutions, insertions, deletions, bulged_start, bulged_end, realigned_target = \
        '', '', '', '', '', '', '', '', 'none'
    if chosen_alignment_b:
        realigned_target, bulged_offtarget_sequence = realignedSequences(targetsite_sequence, chosen_alignment_b, max_score)
        if bulged_offtarget_sequence:
            length = len(chosen_alignment_b.group())
            substitutions, insertions, deletions = chosen_alignment_b.fuzzy_counts
            score = substitutions + (insertions + deletions) * 3
            bulged_start = chosen_alignment_b.start()
            bulged_end = chosen_alignment_b.end()
        else:
            chosen_alignment_strand_b = ''

    return [offtarget_sequence_no_bulge, mismatches, len(offtarget_sequence_no_bulge), chosen_alignment_strand_m, start_no_bulge, end_no_bulge,
            realigned_target,
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end]


""" Get sequences from some reference genome
"""
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][int(start):int(end)]
    elif strand == "-":
        seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
    return str(seq)


def compare(ref, bam, control, targetsite, search_radius, windowsize, mapq_threshold, gap_threshold, start_threshold, mismatch_threshold, name,
            cells, out, all_chromosomes, merged=True):

    output_list = list()

    reference_genome = pyfaidx.Fasta(ref)
    pattern = regex.compile("^[^_]*[0-9XYM]*[^_]*$")

    combined_ga = HTSeq.GenomicArray("auto", stranded=False)  # Store the union of control and nuclease positions
    offtarget_ga_windows = HTSeq.GenomicArray("auto", stranded=False)  # Store potential off-target sites
    ga_narrow_windows = HTSeq.GenomicArray("auto", stranded=False)  # Store potential off-target sites narrow windows read counts

    bg_position = list()  # List to store nuclease_position_counts that were observed at least once
    bg_narrow = list()  # List to store the sum of nuclease_position_counts in the narrow window

    output_folder = os.path.dirname(out)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    output_filename = out + '_count.txt'
    with open(output_filename, 'w') as o:
        print("Writing counts to {0}".format(output_filename), file=sys.stderr)
        if merged:
            print("Tabulate nuclease merged start positions.", file=sys.stderr)
            nuclease_ga, nuclease_ga_windows, nuclease_ga_stranded, nuclease_ga_coverage, total_nuclease_count = \
                tabulate_merged_start_positions(bam, cells, name, targetsite, mapq_threshold, gap_threshold,
                                                start_threshold, out + '_NUCLEASE', pattern, all_chromosomes)
            print("Tabulate control merged start positions.", file=sys.stderr)
            control_ga, control_ga_windows, control_ga_stranded, control_ga_coverage, total_control_count = \
                tabulate_merged_start_positions(control, cells, name, targetsite, mapq_threshold, gap_threshold,
                                                start_threshold, out + '_CONTROL', pattern, all_chromosomes)
        else:
            print("Tabulate nuclease standard start positions.", file=sys.stderr)
            nuclease_ga, nuclease_ga_windows, nuclease_ga_stranded, nuclease_ga_coverage, total_nuclease_count = \
                tabulate_start_positions(bam, cells, name, targetsite, mapq_threshold, gap_threshold, out + '_NUCLEASE', pattern, all_chromosomes)
            print("Tabulate control standard start positions.", file=sys.stderr)
            control_ga, control_ga_windows, control_ga_stranded, control_ga_coverage, total_control_count = \
                tabulate_start_positions(control, cells, name, targetsite, mapq_threshold, gap_threshold, out + '_CONTROL', pattern, all_chromosomes)

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
                    window = HTSeq.GenomicInterval(position.chrom, max(0, position.pos - windowsize),
                                                   position.pos + windowsize + 1)

                    # Start mapping positions, at the specific base position
                    nuclease_position_counts = nuclease_ga[position]
                    control_position_counts = control_ga[position]
                    # Store control_position_counts for which it was observed at least one read
                    if control_position_counts >= 0:
                        bg_position.append(control_position_counts)

                    # In the narrow (parameter-specified) window
                    nuclease_window_counts = sum(nuclease_ga[window])
                    control_window_counts = sum(control_ga[window])

                    # Store control_window_counts greater than zero
                    if control_window_counts >= 0:
                        bg_narrow.append(control_window_counts)

                    # A list of the outputs
                    row = [position.chrom, position.pos, nuclease_position_counts, control_position_counts,
                           nuclease_window_counts, control_window_counts]
                    output_list.append(row)

        print('#Chromosome', 'zero_based_Position', 'Nuclease_Position_Reads', 'Control_Position_Reads',
              'Nuclease_Window_Reads', 'Control_Window_Reads',
              'p_Value', 'narrow_p_Value', 'control_p_Value', 'control_narrow_p_Value', file=o, sep='\t')

        # Empirical cdf
        print(bg_position, bg_narrow)
        ecdf_pos = ECDF(bg_position)
        ecdf_nar = ECDF(bg_narrow)

        # Genomic array to store the p-values for every chromosome:position object
        ga_pval = HTSeq.GenomicArray("auto", typecode='O', stranded=False)

        # Ratio to be used in scaling the nuclease count
        scale_factor = total_control_count/float(total_nuclease_count)

        for idx, fields in enumerate(output_list):
            position_p_val = 1 - ecdf_pos(fields[2]*scale_factor)
            narrow_p_val = 1 - ecdf_nar(fields[4]*scale_factor)

            control_position_p_val = 1 - ecdf_pos(fields[3])
            control_narrow_p_val = 1 - ecdf_nar(fields[5])

            if narrow_p_val < 0.01 or position_p_val < 0.01:
                read_chr = fields[0]
                read_position = fields[1]
                offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
                ga_narrow_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[4]

            print(*(fields + [position_p_val, narrow_p_val, control_position_p_val, control_narrow_p_val]),
                  file=o, sep='\t')

            chr_pos = HTSeq.GenomicPosition(fields[0], int(fields[1]), '.')
            ga_pval[chr_pos] = [position_p_val, narrow_p_val, control_position_p_val, control_narrow_p_val]

        ga_consolidated_windows = find_windows(offtarget_ga_windows, windowsize)    # consolidate windows within 3 bp

        output_alignments(ga_narrow_windows, ga_consolidated_windows, reference_genome, targetsite, name, cells, bam,
                          mismatch_threshold, ga_pval, search_radius, out)

def main():
    parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
    parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
    parser.add_argument('--bam', help='Sorted BAM file', required=True)
    parser.add_argument('--control', help='Control BAM file', required=True)
    parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
    parser.add_argument('--search_radius', help='Search radius around the position window', default=20, type=int)
    parser.add_argument('--windowsize', help='Windowsize', default=3, type=int)
    parser.add_argument('--mapq', help='mapq threshold', default=50, type=int)
    parser.add_argument('--gap', help='Gap threshold', default=3, type=int)
    parser.add_argument('--start', help='Start threshold', default=1 , type=int)
    parser.add_argument('--mismatch_threshold', help='Maximum score threshold', default=6, type=int)
    parser.add_argument('--merged', dest='merged', action='store_true', default=True)
    parser.add_argument('--all_chromosomes', dest='all_chromosomes', action='store_true', default=False)
    parser.add_argument('--name', help='Targetsite Name', required=False)
    parser.add_argument('--cells', help='Cells', required=False)
    parser.add_argument('--out', help='Output file base', required=True)
    args = parser.parse_args()

    # Run the comparison if the control bam is specified, otherwise run the standard site identification routine.
    print("Nuclease: {0}\nControl: {1}".format(args.bam, args.control), file=sys.stderr)
    compare(args.ref, args.bam, args.control, args.targetsite, args.search_radius, args.windowsize, args.mapq, args.gap,
            args.start, args.mismatch_threshold, args.name, args.cells, args.out, args.all_chromosomes, args.merged)

if __name__ == "__main__":
    main()
