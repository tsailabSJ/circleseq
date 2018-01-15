"""
alignReads
"""

from __future__ import print_function

import subprocess
import os
import logging

logger = logging.getLogger('root')
logger.propagate = False

def alignReads(BWA_path, HG19_path, read1, read2, outfile):

    sample_name = os.path.basename(outfile).split('.')[0]
    output_folder = os.path.dirname(outfile)
    base_name = os.path.join(output_folder, sample_name)
    sam_filename = outfile
    bam_filename = base_name + '.bam'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Check if genome is already indexed by bwa
    index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']

    genome_indexed = True
    for extension in index_files_extensions:
        if not os.path.isfile(HG19_path + extension):
            genome_indexed = False
            break

    # If the genome is not already indexed, index it
    if not genome_indexed:
        logger.info('Genome index files not detected. Running BWA to generate indices.')
        bwa_index_command = '{0} index {1}'.format(BWA_path, HG19_path)
        logger.info('Running bwa command: %s', bwa_index_command)
        subprocess.call(bwa_index_command.split())
        logger.info('BWA genome index generated')
    else:
        logger.info('BWA genome index found.')

    # Run paired end alignment against the genome
    logger.info('Running paired end mapping for {0}'.format(sample_name))
    bwa_alignment_command = '{0} mem {1} {2} {3} > {4}'.format(BWA_path, HG19_path, read1, read2, sam_filename)
    samtools_sam_to_bam_command = 'samtools sort -o {0} {1}'.format(bam_filename, sam_filename)
    samtools_index_command = 'samtools index {0}'.format(bam_filename)
    samtools_sort_by_name_command = 'samtools sort -o {0} -n {1}'.format("".join([base_name, '_sorted.bam']), bam_filename)

    # Open the outfile and redirect the output of the alignment to it.
    logger.info(bwa_alignment_command)
    subprocess.check_call(bwa_alignment_command, shell=True)
    logger.info('Paired end mapping for {0} completed.'.format(sample_name))

    # Convert SAM to BAM file
    logger.info(samtools_sam_to_bam_command)
    subprocess.check_call(samtools_sam_to_bam_command, shell=True)
    logger.info('Sorting by coordinate position for {0} complete.'.format(sample_name))

    # Index BAM file
    logger.info(samtools_index_command)
    subprocess.check_call(samtools_index_command, shell=True)
    logger.info('Indexing for {0} complete.'.format(sample_name))

    # Sort BAM file by name
    logger.info(samtools_sort_by_name_command)
    subprocess.check_call(samtools_sort_by_name_command, shell=True)
    logger.info('Sorting for {0} by name complete.'.format(sample_name))




