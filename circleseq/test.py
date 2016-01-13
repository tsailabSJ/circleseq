import regex

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
        pattern = '(?V1)(?b:' + pattern + ')'
    if mismatches > 0:
        pattern = pattern + '{{i<=1,d<=1,s<=6,e<={0}}}'.format(mismatches)
    return pattern

"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match.
"""
def alignSequences(targetsite_sequence, window_sequence, max_mismatches=6):
    # Try both strands
    query_regex = regexFromSequence(targetsite_sequence, mismatches=max_mismatches)
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

        return [match_sequence, distance, length, strand, start, end]

def main():
    target = 'GAGTCCGAGCAGAAGAAGAANGG'
    # windowsequence = 'GAAGGGCCTGAGTCCGAGCAGAAGAAGAAGGGCTCCCATCACATCAACCGGT'
    windowsequence = 'GGCCTGAGTCCGAGCAGAAGAAGAAGGGCTCCCATCACATCAAC'

    pattern = regexFromSequence(target, mismatches=7)

    match = regex.search(pattern, windowsequence, flags=regex.BESTMATCH)

    matches = regex.findall(pattern, windowsequence)

    print(matches)
if __name__ == "__main__":
    main()