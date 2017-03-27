#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_circleseq
----------------------------------

Tests for `circleseq` module.
"""

import yaml
import unittest
import os
import shutil
import utils
from circleseq import circleseq

TEST_SAMPLE_NAME = 'EMX1'
TEST_OUTPUT_PATH = 'test_output'
TEST_MIN_READS = 1000
TEST_DEMULTIPLEX_MANIFEST_PATH = os.path.join(TEST_OUTPUT_PATH, 'demultiplex_manifest.yaml')
TEST_MANIFEST_PATH = os.path.join(TEST_OUTPUT_PATH, 'test_manifest.yaml')

TEST_BWA_PATH = 'bwa'
TEST_BEDTOOLS_PATH = 'bedtools'

TEST_REFERENCE_GENOME = 'test_genome.fa'

CORRECT_ALIGNED_OUTPUT = 'data/aligned'
CORRECT_IDENTIFIED_OUTPUT = 'data/identified'
CORRECT_MERGED_OUTPUT = 'data/merged'
CORRECT_VISUALIZATION_OUTPUT = 'data/visualization'

CORRECT_ALL_OUTPUT = 'data'

class FullPipelineTestCase(unittest.TestCase):

    def setUp(self):
        # Create the test output folder
        os.makedirs(TEST_OUTPUT_PATH)

        # Create the test demultiplexing YAML
        test_manifest_data = {}
        test_manifest_data['undemultiplexed'] = TEST_UNDEMULTIPLEXED_FILES
        test_manifest_data['demultiplex_min_reads'] = TEST_MIN_READS
        test_manifest_data['samples'] = TEST_SAMPLES
        test_manifest_data['output_folder'] = TEST_OUTPUT_PATH
        test_manifest_data['bwa'] = TEST_BWA_PATH
        test_manifest_data['bedtools'] = TEST_BEDTOOLS_PATH
        test_manifest_data['reference_genome'] = TEST_REFERENCE_GENOME

        with open(TEST_MANIFEST_PATH, 'w') as f:
            f.write(yaml.dump(test_manifest_data, default_flow_style=False))


    def testFullPipeline(self):
        c = circleseq.CircleSeq()
        c.parseManifest(TEST_MANIFEST_PATH)

        # Align and test the alignment output
        c.alignReads()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'aligned'), CORRECT_ALIGNED_OUTPUT))

        # Identify offtargets and test the output
        c.identifyOfftargetSites()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'identified'), CORRECT_IDENTIFIED_OUTPUT))

        # Merge identified sites
        c.mergeAlignReads()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'merged'), CORRECT_MERGED_OUTPUT))

        # Visualize filtered sites
        c.visualize()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'visualization'), CORRECT_VISUALIZATION_OUTPUT))


    def tearDown(self):
        # Delete temp output
        shutil.rmtree(TEST_OUTPUT_PATH)
        pass

if __name__ == '__main__':
    unittest.main()