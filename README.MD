CIRCLE-seq: Comprehensive In vitro Reporting of CLeavage Effects by Sequencing
==================================================================================

This is a repository for CIRCLE-seq analytical software, which takes sample-specific paired end FASTQ files as input and produces a list of CIRCLE-seq detected off-target cleavage sites as output.

Dependencies
=======

* Python (2.7)
* [bwa] alignment tool
* [samtools] alignment format manipulation

[bwa]: http://bio-bwa.sourceforge.net
[samtools]: http://www.htslib.org

**Usage**
=========

Using this tool is simple, just create a ``.yaml`` manifest file referencing the dependencies and sample ``.fastq.gz`` file paths. Then, run ``python /path/to/circleseq.py all --manifest /path/to/manifest.yaml``

Below is an example ``manifest.yaml`` file::

    reference_genome: /data/joung/genomes/Homo_sapiens_assembly19.fasta
    analysis_folder: /data/joung/CIRCLE-Seq/test2
    
    bwa: bwa
    samtools: samtools
    
    samples:
      VEGFA_site1:
          target: GGGTGGGGGGAGTTTGCTCCNGG
          read1: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/1_S1_L001_R1_001.fastq.gz
          read2: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/1_S1_L001_R2_001.fastq.gz
          description: VEGFA_site1
      EMX1:
          target: GAGTCCGAGCAGAAGAAGAANGG
          read1: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/2_S2_L001_R1_001.fastq.gz
          read2: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/2_S2_L001_R2_001.fastq.gz
          description: EMX1

