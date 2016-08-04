import string
import re
import gzip
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
            l1 = f.readline().rstrip('\n')
            if not l1:
                break
            l2 = f.readline().rstrip('\n')
            l3 = f.readline().rstrip('\n')
            l4 = f.readline().rstrip('\n')
            yield [l1, l2, l3, l4]

def reverseComplement(sequence):
    transtab = string.maketrans("ACGT","TGCA")
    return sequence.translate(transtab)[::-1]
