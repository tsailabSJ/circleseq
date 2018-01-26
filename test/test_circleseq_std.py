#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_circleseq_std
----------------------------------

Tests for `circleseq` module.
"""

import yaml
import unittest
import os
import shutil
import utils
from circleseq import circleseq

TEST_OUTPUT_PATH = 'tmp'

TEST_MANIFEST_PATH = os.path.join('CIRCLEseq_StandardTest.yaml')

CORRECT_ALIGNED_OUTPUT = 'data/StandardOutput/aligned'
CORRECT_IDENTIFIED_OUTPUT = 'data/StandardOutput/identified'
CORRECT_VARIANTS_OUTPUT = 'data/StandardOutput/variants'
CORRECT_VISUALIZATION_OUTPUT = 'data/StandardOutput/visualization'

CORRECT_ALL_OUTPUT = 'data'

class FullPipelineTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def testFullPipeline(self):
        c = circleseq.CircleSeq()
        c.parseManifest(TEST_MANIFEST_PATH)
        
        # Align and test the alignment output
        c.alignReads()
        self.assertTrue(utils.checkFolderEquality(os.path.join(c.analysis_folder, "aligned"), CORRECT_ALIGNED_OUTPUT))

        # Find cleavage sites
        c.findCleavageSites()
        self.assertTrue(utils.checkFolderEquality(os.path.join(c.analysis_folder, 'identified'), CORRECT_IDENTIFIED_OUTPUT))

        # Visualize filtered sites
        c.visualize()
        self.assertTrue(utils.checkFolderEquality(os.path.join(c.analysis_folder, 'visualization'), CORRECT_VISUALIZATION_OUTPUT))

        # Look for genomic variants
        c.callVariants()
        self.assertTrue(utils.checkFolderEquality(os.path.join(c.analysis_folder, 'variants'), CORRECT_VARIANTS_OUTPUT))


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()