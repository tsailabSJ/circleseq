# CIRCLE-seq: Circularization In vitro Reporting of CLeavage Effects by Sequencing

This is a repository for CIRCLE-seq analytical software, which takes sample-specific paired end FASTQ files as input and produces a list of CIRCLE-seq detected off-target cleavage sites as output.

## Table of Contents
- [Features](#features)
- [Dependencies](#dependencies)
- [Getting Set Up](#setup)
	- [Installing Dependencies](#install_dependencies)
	- [Download Reference Genome](#reference_genome)
	- [Download and Set Up circleseq](#circleseq_setup)
- [Usage](#usage)
	- [Quickstart](#quickstart)
	- [Writing A Manifest File](#write_manifest)
	- [Pipeline Outputs](#pipeline_output)
- [Testing the circleseq Package](#testing)
- [Frequently Asked Questions](#FAQ)

## Features<a name="features"></a>

This package implements a pipeline that takes in reads from the CIRCLE-seq assay and returns detected cleavage sites as output. The individual pipeline steps are:

1. **Merge**: Merge read1 an read2 for easier mapping to genome.
2. **Read Alignment**: Merged paired end reads from the assay are aligned to the reference genome using the BWA-MEM algorithm with default parameters (Li. H, 2009).
3. **Cleavage Site Identification**: Mapped sites are analyzed to determine which represent high-quality cleavage sites.
4. **Visualization of Results**: Identified on-target and off-target cleavage sites are rendered as a color-coded alignment map for easy analysis of results.

## Dependencies<a href="dependencies"></a>

* Python (2.7)
* Reference genome fasta file ([Example](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta))
* [`bwa`](https://github.com/lh3/bwa) alignment tool
* [`samtools`](http://www.htslib.org/download/) alignment format manipulation


## Getting Set Up<a href="setup"></a>


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

1. Download a copy of the circleseq package source code. You can either download and unzip the latest source from the github [release page](https://github.com/tsailabSJ/circleseq/releases), or you use git to clone the repository by running `git clone --recursive https://github.com/tsailabSJ/circleseq.git`
2. Install circleseq dependencies by entering the circleseq directory and running `pip install -r requirements.txt`

Once all circleseq dependencies are installed, you will be ready to start using circleseq.


## Usage<a name="usage"></a>

### Quickstart<a name="quickstart"></a>
Using this tool is simple. After [getting set up](#setup), create a ``.yaml`` manifest file referencing the dependencies and sample ``.fastq.gz`` file paths. Then, run ``python /path/to/circleseq.py all --manifest /path/to/manifest.yaml``

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


### Writing A Manifest File<a name="write_manifest"></a>
When running the end-to-end analysis functionality of the circleseq package a number of inputs are required. To simplify the formatting of these inputs and to encourage reproducibility, these parameters are inputted into the pipeline via a manifest formatted as a YAML file. YAML files allow easy-to-read specification of key-value pairs. This allows us to easily specify our parameters. The following fields are required in the manifest:

- `reference_genome`: The absolute path to the reference genome FASTA file.
- `output_folder`: The absolute path to the folder in which all pipeline outputs will be saved.
- `bwa`: The absolute path to the `bwa` executable
- `samtools`: The absolute path to the `samtools` executable
- `read_threshold`: The minimum number of reads at a location for that location to be called as a site. We recommend leaving it to the default value of 4.
- `window_size`: Size of the sliding window, we recommend leaving it to the default value of 3.
- `mapq_threshold`: Minimum read mapping quality score. We recommend leaving it to the default value of 50.
- `start_threshold`: Tolerance for breakpoint location. We recommend leaving it to the default value of 1.
- `gap_threshold`: Number of tolerated gaps in the fuzzy target search setp. We recommend leaving it to the default value of 3.
- `mismatch_threshold`: Number of tolerated gaps in the fuzzy target search setp. We recommend leaving it to the default value of 6.
- `merged_analysis`: Whether or not the paired read merging step should takingTrue
- `samples`: Lists the samples you wish to analyze and the details for each. Each sample name should be nested under the top level samples key, and each sample detail should be nested under the sample name. See the sample manifest for an example.
    - For each sample, you must provide the following parameters:
        - `target`: Target sequence for that sample. Accepts degenerate bases.
        - `read1`: The absolute path to the .FASTQ(.gz) file containing the read1 reads.
        - `read2`: The absolute path to the .FASTQ(.gz) file containing the read2 reads.
        - `controlread1`: The absolute path to the .FASTQ(.gz) file containing the control read1 reads.
        - `controlread2`: The absolute path to the .FASTQ(.gz) file containing the control read2 reads.
        - `description`: A brief description of the sample


### Pipeline Output<a name="pipeline_output"></a>
When running the full pipeline, the results of each step are outputted to the `output_folder` in a separate folder for each step. The output folders and their respective contents are as follows:

- `output_folder/aligned`: Contains an alignment `.sam`, alignment `.bam`, sorted `bam`, and `.bai` index file for each sample.
- `output_folder/fastq`: Merged `.fastq.gz` files for each sample.
- `output_folder/identified`: Contains tab-delimited `.txt` files for each sample containing the identified DSBs, control DSBs, filtered DSBs, and read quantification.
- `output_folder/visualization`: Contains a `.svg` vector image representing an alignment of all detected off-targets to the targetsite for each sample.


## Testing the circleseq Package<a name="testing"></a>

In the spirit of Test-Driven Development, we have written tests for the pipeline to protect against regressions when making future updates. These can be used to ensure that the software is running with expected functionality.

NOTE: Due to differences in sorting between different versions of the `bwa` package, you must be using `bwa 0.7.11` for these tests to work. We also recommend that you use `samtools 1.3` when running these tests for consistency's sake.

### Single-Step Regression Tests<a name="regression_tests"></a>

For ongoing testing and development, we have created an abridged set of input data and expected output data for each step of the pipeline. This way, changes to the pipeline can be quickly tested for feature regression.

To run these tests, you must first install the `nose` testing Python package.

```
pip install nose
```

Then, from the circleseq root directory, simply run

```
nosetests
```

and the regression tests for for the pipeline will be run.


## FAQ<a name="FAQ"></a>

None yet, we will keep this updated as needed.
