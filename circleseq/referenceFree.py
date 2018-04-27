from __future__ import print_function

import argparse
import itertools
import regex
import re
import gzip
import sys
import collections
from findCleavageSites import regexFromSequence, alignSequences, reverseComplement, extendedPattern, realignedSequences

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
Main function to find off-target sites in reference-free fashion
"""
def analyze(fastq1_filename, fastq2_filename, targetsite, out_base, name='', cells='', mismatch_threshold=7):

    read_count = 0
    c = collections.Counter()
    d = collections.defaultdict(list)

    fastq1_file = fq(fastq1_filename)
    fastq2_file = fq(fastq2_filename)
    for r1, r2 in itertools.izip(fastq1_file, fastq2_file):
        r1_sequence = r1[1].rstrip('\n')
        r2_sequence = r2[1].rstrip('\n')
        joined_seq = reverseComplement(r1_sequence) + r2_sequence
        truncated_joined_seq = joined_seq[130:170]

        sequence_data = alignSequences(targetsite, truncated_joined_seq, max_score=mismatch_threshold)
        offtarget, mismatch, length, strand, start, end, realigned_target = sequence_data[:7]

        if offtarget:
            c[offtarget] += 1
            d[offtarget].append(joined_seq)

        read_count += 1
        if not read_count % 100000:
            print(read_count/float(1000000), end=" ", file=sys.stderr)

    print('Finished tabulating reference-free discovery counts.', file=sys.stderr)
    out_filename = out_base + '.txt'

    with open(out_filename, 'w') as o:
        for target_sequence, target_count in c.most_common():
            print(target_sequence, target_count, file=o)
            off_target_fasta_filename = '{0}_{1:04d}_{2}.fasta'.format(out_base, target_count, target_sequence)
            with open(off_target_fasta_filename, 'w') as off_target_fasta_file:
                j = 0
                for sequence in d[target_sequence]:
                    j += 1
                    print('>{0:04d}_{1}_{2}'.format(target_count, target_sequence, j), file=off_target_fasta_file)
                    print(sequence, file=off_target_fasta_file)

def join_write_output(fastq1_filename, fastq2_filename, out):
    fastq1_file = fq(fastq1_filename)
    fastq2_file = fq(fastq2_filename)

    with open(out, 'w') as o:
        for r1, r2 in itertools.izip(fastq1_file, fastq2_file):
            header = '>{0}'.format(r1[0])
            r1_sequence = r1[1].rstrip('\n')
            r2_sequence = r2[1].rstrip('\n')
            joined_seq = reverseComplement(r1_sequence) + r2_sequence
            print(header, end='', file=o)
            print(joined_seq, file=o)


def main():
    parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
    parser.add_argument('--fq1', help='FASTQ Read 1', required=True)
    parser.add_argument('--fq2', help='FASTQ Read 2', required=True)
    parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
    parser.add_argument('--name', help='Targetsite Name', required=False)
    parser.add_argument('--cells', help='Cells', required=False)
    parser.add_argument('--mismatch_threshold', help='Maximum score threshold', default=7, type=int)
    parser.add_argument('--out', help='Output file base', required=True)
    args = parser.parse_args()

    analyze(args.fq1, args.fq2, args.targetsite, args.out, args.name, args.cells, args.mismatch_threshold)

if __name__ == "__main__":
    main()
