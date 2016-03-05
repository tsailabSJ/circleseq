import regex
import nwalign as nw
import swalign
import string

def reverseComplement(sequence):
    transtab = string.maketrans("ACGT","TGCA")
    return sequence.translate(transtab)[::-1]

def regexFromSequence(seq, lookahead=True, indels=1, errors=7):
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
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match.
"""
def alignSequences(targetsite_sequence, window_sequence, max_mismatches=7):
    # Try both strands
    query_regex_standard, query_regex_gap = regexFromSequence(targetsite_sequence, errors=max_mismatches)

    alignments = list()
    alignments.append(('+', 'standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
    alignments.append(('-', 'standard', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
    alignments.append(('+', 'gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
    alignments.append(('-', 'gapped', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))

    top_distance_score = 0
    chosen_alignment = None
    for i, aln in enumerate(alignments):
        strand, alignment_type, match = aln
        if match != None:
            substitutions, insertions, deletions = match.fuzzy_counts
            distance_score = substitutions + (insertions + deletions) * 3
            if distance_score > top_distance_score:
                chosen_alignment = match
                top_distance_score = distance_score
                print(match, distance_score)

    if chosen_alignment:
        match_sequence = chosen_alignment.group()
        distance = sum(chosen_alignment.fuzzy_counts)
        length = len(match_sequence)
        start = chosen_alignment.start()
        end = chosen_alignment.end()
        return [match_sequence, distance, length, strand, start, end]
    else:
        return [''] * 6




    # if forward_alignment is None and reverse_alignment is None:
    #     return ['', '', '', '', '', '']
    # else:
    #     if forward_alignment is None and reverse_alignment is not None:
    #         strand = '-'
    #         alignment = reverse_alignment
    #     elif reverse_alignment is None and forward_alignment is not None:
    #         strand = '+'
    #         alignment = forward_alignment
    #     elif forward_alignment is not None and reverse_alignment is not None:
    #         forward_distance = sum(forward_alignment.fuzzy_counts)
    #         reverse_distance = sum(reverse_alignment.fuzzy_counts)
    #
    #         if forward_distance > reverse_distance:
    #             strand = '-'
    #             alignment = reverse_alignment
    #         else:
    #             strand = '+'
    #             alignment = forward_alignment
    #
    #     match_sequence = alignment.group()
    #     distance = sum(alignment.fuzzy_counts)
    #     length = len(match_sequence)
    #     start = alignment.start()
    #     end = alignment.end()
    #
    #     return [match_sequence, distance, length, strand, start, end]

def alignSequences2(ref_seq, query_seq):
    match = 2
    mismatch = -1
    ref_length = len(ref_seq)
    matches_required = len(ref_seq) - 1 - 7  # allow up to 8 mismatches
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=-3, gap_extension_penalty=-100, prefer_gap_runs=True)  # you can also choose gap penalties, etc...
    # sw = swalign.LocalAlignment(scoring, gap_penalty=-10, gap_extension_penalty=-0.5, prefer_gap_runs=True)  # you can also choose gap penalties, etc...
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
    # target = 'TTTNCTGATGGTCCATGTCTGTTACTC'


    # windowsequence = 'AATGTGTGTCTGCTGGAAGCTCCTATTCTTCCGCCATTTTCCAGTCCTCCAGAAGTTTCCTGATGGTCCATGTCTGAATTAGACACCCCTCTTCTTTGTTCCAGTTGCACCTGTAATTCTTCAGCATAGTACTTCTTAAACTGTTTTTAA'
    # windowsequence = 'GGCCTGAGTCCGAGCAGAAGCAAGAAGGGCTCCCATCACATCAAC'

    target = 'TTTNGGGACGGGGAGAAGGAAAAGAGG'
    windowsequence = 'AATTTGGGGGGATTCATTACTCTATTTGGATTTGTTAGGGAGGAAGGCAGGTGGGATTTTTCTTCTCATTCTTATCTCTTTCCTTCTTCCCGTCCCAGAAAGAAACTAAGAATAATAACCAAATTATTAAAATGACTCACCGCCCTTCCA'

    print(alignSequences(target, windowsequence, max_mismatches=7))


if __name__ == "__main__":
    main()