import nwalign as nw
import Levenshtein as l
import difflib
import os

def main():

    # a = 'GCCTGAGTCCGAGCAGAAGAAGAAGGGCTCCCATCACATCAAC'
    # b = 'GAGTCGAGCAGAAGAAGAANGG'

    a = 'AATGTGTGTCTGCTGGAAGCTCCTATTCTTCCGCCATTTTCCAGTCCTCCAGAAGTTTCCTGATGGTCCATGTCTGAATTAGACACCCCTCTTCTTTGTTCCAGTTGCACCTGTAATTCTTCAGCATAGTACTTCTTAAACTGTTTTTAA'
    b= 'TTTNCTGATGGTCCATGTCTGTTACTC'

    print(l.distance(a, b))
    print(l.editops(a, b))
    print(l.matching_blocks(l.editops(a,b), a, b))



if __name__ == "__main__":
    main()