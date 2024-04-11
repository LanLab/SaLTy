from time import sleep as sl
import sys
import unittest
from salty import salty
import os
import argparse
import shutil

class TestSalty_full(unittest.TestCase):
    def setUp(self):
        cwd = os.getcwd()
        self.outputfolder = cwd+"/salty/tests/outputs"
        if os.path.exists(self.outputfolder):
            shutil.rmtree(self.outputfolder)
            os.mkdir(self.outputfolder)
            sl(0.5)
        else:
            os.mkdir(self.outputfolder)
            sl(0.5)

    def test_with_genomes(self):

        self.maxDiff = 5000
        cwd = os.getcwd()
        args = argparse.Namespace(input_folder = cwd+"/dodge/tests/inputs/",
                                    threads = 1,
                                    force = False,
                                    summary = True,
                                    report=False,
                                    version=False,
                                    check=False,
                                    output_folder=self.outputfolder+"/",
                                    csv_format = False,
                                    lineages = cwd+"/salty/resources/alleles/alleles.csv",
                                    kma_index = cwd+"/salty/resources/kmaIndex/kmaIndex",
                                    mlstPrediction = True
                                    )

        salty.main(args)

        outsummary = self.outputfolder+"/summaryReport.csv"

        expected_outsummary = cwd+"/salty/tests/expected_outputs/summaryReport.csv"


        self.assertListEqual(
            list(open(outsummary)),
            list(open(expected_outsummary)))

    def tearDown(self):
        shutil.rmtree(self.outputfolder)