from __future__ import print_function

import argparse

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
def analyze(fastq1_filename, fastq2_filename, targetsite, name, cells, out):
    fastq1_file = fq(fastq1_filename)
    fastq2_file = fq(fastq2_filename)
    pass



def main():
    parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
    parser.add_argument('--fq1', help='FASTQ Read 1', required=True)
    parser.add_argument('--fq1', help='FASTQ Read 2', required=True)
    parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
    parser.add_argument('--name', help='Targetsite Name', required=False)
    parser.add_argument('--cells', help='Cells', required=False)
    parser.add_argument('--out', help='Output file base', required=False)
    args = parser.parse_args()

    analyze(args.bam, args.targetsite, args.name, args.cells, args.out)

if __name__ == "__main__":
    main()
