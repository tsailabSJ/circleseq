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
    
    read_threshold: 4
    window_size: 3
    mapq_threshold: 50
    start_threshold: 1
    gap_threshold: 3
    mismatch_threshold: 6
    merged_analysis: True
    
    samples:
        U2OS_exp1_VEGFA_site_1:
            target: GGGTGGGGGGAGTTTGCTCCNGG
            read1: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/1_S1_L001_R1_001.fastq.gz
            read2: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/1_S1_L001_R2_001.fastq.gz
            controlread1: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/4_S4_L001_R1_001.fastq.gz
            controlread2: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/4_S4_L001_R2_001.fastq.gz
            description: U2OS_exp1
        U2OS_exp1_EMX1:
            target: GAGTCCGAGCAGAAGAAGAANGG
            read1: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/2_S2_L001_R1_001.fastq.gz
            read2: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/2_S2_L001_R2_001.fastq.gz
            controlread1: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/4_S4_L001_R1_001.fastq.gz
            controlread2: /data/joung/sequencing_fastq/150902_M01326_0235_000000000-AHLT8/fastq/4_S4_L001_R2_001.fastq.gz
            description: U2OS_exp1

**Getting Set Up**
=========

### Install Dependencies<a name="install_dependencies"></a>

To run circleseq, you must first install all necessary dependencies:

- **Python 2.7**: If a version does not come bundled with your operating system, we recommend the [Anaconda](https://www.continuum.io/downloads) scientific Python package.
- **Burrows-Wheeler Aligner (bwa)**: You can either install bwa with a package manager (e.g. `brew` on OSX or `apt-get` on Ubuntu/Debian), or you can download it from the [project page](http://bio-bwa.sourceforge.net/) and compile it from source.
- **Samtools**: You can either install samtools with a package manager (e.g. `brew` or `apt-get`), or you can download it from the [project page](http://www.htslib.org/download/) and compile it from source.

For both bwa and samtools, make sure you know the path to the respective executables, as they need to be specified in the pipeline manifest file.

### Download Reference Genome<a name="reference_genome"></a>

The circleseq package requires a reference genome for read mapping. You can use any genome of your choosing, but for all of our testing and original CIRCLE-seq analyses we use hg19 ([download](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta)). Be sure to (g)unzip the FASTA file before use if it is compressed.

### Download and Set Up circleseq<a name="circleseq_setup"></a>

Once all dependencies are installed, there are a few easy steps to download and set up the circleseq package:

1. Obtain a copy of the circleseq package source code. You can either download and unzip the latest source from the github [release page](https://github.com/tsailabSJ/circleseq/releases), or you use git to clone the repository by running `git clone --recursive https://github.com/tsailabSJ/circleseq.git`
2. Install guideseq dependencies by entering the circleseq directory and running `pip install -r requirements.txt`

Once all guideseq dependencies are installed, you will be ready to start using circleseq.