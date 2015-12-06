"""
circleseq.py as the wrapper for CIRCLE-seq analysis
"""

from alignReads import alignReads

import argparse
import os
import sys
import traceback
import log
import yaml
import validation
import findCleavageSites

logger = log.createCustomLogger('root')

class CircleSeq:

    def __init__(self):
        pass

    def parseManifest(self, manifest_path):
        logger.info('Loading manifest...')

        with open(manifest_path, 'r') as f:
            manifest_data = yaml.load(f)

        try:
            # Validate manifest data
            validation.validateManifest(manifest_data)

            self.BWA_path  = manifest_data['bwa']
            self.reference_genome = manifest_data['reference_genome']
            self.analysis_folder = manifest_data['analysis_folder']
            self.samples = manifest_data['samples']

        except Exception as e:
            logger.error('Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
            sys.exit()

    def alignReads(self):
        logger.info('Aligning reads...')

        try:
            self.aligned = {}
            self.aligned_sorted = {}
            for sample in self.samples:
                sample_alignment_path = os.path.join(self.analysis_folder, 'aligned', sample + '.sam')
                alignReads(self.BWA_path,
                           self.reference_genome,
                           self.samples[sample]['read1'],
                           self.samples[sample]['read2'],
                           sample_alignment_path)
                self.aligned[sample] = sample_alignment_path
                self.aligned_sorted[sample] = os.path.join(self.analysis_folder, 'aligned', sample + '_sorted.bam')
                logger.info('Finished aligning reads to genome.')

        except Exception as e:
            logger.error('Error aligning')
            logger.error(traceback.format_exc())
            quit()

    def findCleavageSites(self):
        logger.info('Identifying off-target cleavage sites.')
        try:

            for sample in self.samples:
                sorted_bam_file = os.path.join(self.analysis_folder, 'aligned', sample + '_sorted.bam')
                identified_sites_file = os.path.join(self.analysis_folder, 'identified', sample + '_identified.txt')
                findCleavageSites.analyze(self.reference_genome, sorted_bam_file, self.samples[sample]['target'], 4, 3, True, sample, self.samples[sample]['description'], identified_sites_file)
        except Exception as e:
            logger.error('Error identifying off-target cleavage site.')
            logger.error(traceback.format_exc())
            quit()

    def visualize(self):
        pass



def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')

    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)

    align_parser = subparsers.add_parser('align', help='Run alignment only')
    align_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)

    identify_parser = subparsers.add_parser('identify', help='Run identification only')
    identify_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)

    visualize_parser = subparsers.add_parser('identify', help='Run visualization only')
    visualize_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)

    return parser.parse_args()

def main():
    args = parse_args()

    if args.command == 'all':
        c = CircleSeq()
        c.parseManifest(args.manifest)
        c.alignReads()
        c.findCleavageSites()
    elif args.command == 'align':
        c = CircleSeq()
        c.parseManifest(args.manifest)
        c.alignReads()
    elif args.command == 'identify':
        c = CircleSeq()
        c.parseManifest(args.manifest)
        c.findCleavageSites()
    elif args.command == 'visualize':
        c = CircleSeq()
        c.parseManifest(args.manifest)
        c.visualize()

if __name__ == '__main__':
    main()
