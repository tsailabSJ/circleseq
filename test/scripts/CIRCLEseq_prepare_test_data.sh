#####################################################################################
### CIRCLEseq_prepare_test_data.sh
#####################################################################################

### To be run from  '/data/joung/users/jose/CIRCLEseq/FilesTest/aligned'

# Regions to consider
on_target="2:73160981-73161004"
off_target="8:120587494-120587517 11:30301802-30301825"
hotspots="2:73161104-73161159"


# Get the names of reads that overlap with the selected test regions:
samtools view U2OS_exp1_EMX1.bam $on_target $off_target $hotspots | cut -f1 | sort | uniq > read_names.txt


# Subset FASTQs to extract _all_ read pairs where at least one of the reads falls in a specified test region
zcat /data/joung/users/jose/CIRCLEseq/FilesTest/fastq/2_S2_L001_R1_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > /data/joung/users/jose/CIRCLEseq/FilesTest/test/input_data/EMX1.r1.fastq.gz
zcat /data/joung/users/jose/CIRCLEseq/FilesTest/fastq/2_S2_L001_R2_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > /data/joung/users/jose/CIRCLEseq/FilesTest/test/input_data/EMX1.r2.fastq.gz
zcat /data/joung/users/jose/CIRCLEseq/FilesTest/fastq/4_S4_L001_R1_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > /data/joung/users/jose/CIRCLEseq/FilesTest/test/input_data/EMX1_control.r1.fastq.gz
zcat /data/joung/users/jose/CIRCLEseq/FilesTest/fastq/4_S4_L001_R2_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > /data/joung/users/jose/CIRCLEseq/FilesTest/test/input_data/EMX1_control.r2.fastq.gz
