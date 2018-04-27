#!/usr/bin/env bash
#####################################################################################
### CIRCLEseq_prepare_test_genome.sh: assemble reference test
#####################################################################################
### Get chromosomes
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.1.fa.gz
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.2.fa.gz
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.4.fa.gz
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.8.fa.gz
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.12.fa.gz

### Asemble reference
cat *.fa.gz > Homo_sapiens.GRCh37.subset.fa.gz
gunzip Homo_sapiens.GRCh37.subset.fa.gz
samtools faidx Homo_sapiens.GRCh37.subset.fa

### Pad test regions with 10kb on either side
bedtools slop -i CIRCLEseq_test.bed -g Homo_sapiens.GRCh37.subset.fa.fai -b 10000 > CIRCLEseq_test_padded.bed

### Extract sequences from reference file for each paded interval
bedtools getfasta -fi Homo_sapiens.GRCh37.subset.fa -bed CIRCLEseq_test_padded.bed -fo CIRCLEseq_test_genome.fa -name
