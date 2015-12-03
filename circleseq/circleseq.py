"""
circleseq.py as the wrapper for CIRCLE-seq analysis
"""

import argparse
import log
import sys
import yaml

logger = log.createCustomLogger('root')

class CircleSeq:

    def __init__(self):
        pass


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
            self.bedtools = manifest_data['bedtools']
            self.reference_genome = manifest_data['reference_genome']
            self.output_folder = manifest_data['output_folder']
            self.samples = manifest_data['samples']

        except Exception as e:
            logger.error('Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
            sys.exit()


def main():
    pass

if __name__ == '__main__':
    main()
