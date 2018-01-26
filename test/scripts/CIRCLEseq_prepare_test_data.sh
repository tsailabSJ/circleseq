#!/usr/bin/env bash
#####################################################################################
### CIRCLEseq_prepare_test_data.sh: assemble the fastq files for the test
#####################################################################################
### Regions
on_target="2:73160981-73161004"
off_target01="8:120587494-120587517"
off_target02="1:234492858-234492881"
off_target03="12:73504668-73504691"
off_target04="4:48639390-48639413"
hotspots="1:121485221-121485228"

### Get the names of reads that overlap with the selected test regionsq
samtools view sample.bam $on_target $off_target01 $off_target02 $off_target03 $off_target04 $hotspots | cut -f1 | sort | uniq > sample_read_names.txt
samtools view control.bam $on_target $off_target01 $off_target02 $off_target03 $off_target04 $hotspots | cut -f1 | sort | uniq > control_read_names.txt
cat sample_read_names.txt control_read_names.txt > read_names.txt

### Subset FASTQs to extract _all_ read pairs where at least one of the reads falls in a specified test region
zcat fastq/128_S3_L001_R1_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > TEST.r1.fastq.gz
zcat fastq/128_S3_L001_R2_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > TEST.r2.fastq.gz
zcat fastq/Negative_S1_L001_R1_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > TEST_control.r1.fastq.gz
zcat fastq/Negative_S1_L001_R2_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > TEST_control.r2.fastq.gz
