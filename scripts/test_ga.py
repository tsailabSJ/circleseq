import HTSeq

def main():
    ga = HTSeq.GenomicArray("auto", typecode='O', stranded=False)
    position = HTSeq.GenomicPosition('chr1', 123203, '.')

    ga[HTSeq.GenomicInterval( "chr1", 100000, 101000 , "." )] = [0.05, 0.002, 0.04, 0.005]

    iv = HTSeq.GenomicInterval( "chr1", 100000, 130000 , "." )

    for interval, value in ga[iv].steps():
        print(interval, value)

if __name__ == "__main__":
    main()