
Comprehensive In vitro Reporting of CLeavage Effects by Sequencing (CIRCLE-seq)
===============================================================================
This is a repository for CIRCLE-seq analytical software, which takes sample-specific paired end FASTQ files as input and produces a list of CIRCLE-seq detected off-target cleavage sites as output.

Usage
=====

1. The first step is to map the reads to the human reference genome using bwa.

```
bwa mem /data/joung/genomes/Homo_sapiens_assembly19.fasta 1025_S1_L001_R1_001.fastq.gz 1025_S1_L001_R2_001.fastq.gz > 1025_S1_L001.sam
samtools view -bS 1025_S1_L001.sam | samtools sort - 1025_S1_L001
samtools index 1025_S1_L001.bam
```
2. The second step is to sort the reads by read name to group paired end reads.

```
samtools sort -n 1025_S1_L001.bam 1025_S1_L001_sorted
```

3. The final step is to run FindCleavageSites-Circular.py to identify CIRCLE-seq Cas9 off-target cleavage sites.

```
python /data/joung/bitbucket/circle-seq/FindCleavageSites-Circular.py --ref /data/joung/genomes/Homo_sapiens_assembly19.fasta --reads 4 --windowsize 3 --bam 1025_S1_L001_sorted.bam --targetsite GGGAAAGACCCAGCATCCGTNGG > CIRCLE-Seq_Adli_site1.txt
```