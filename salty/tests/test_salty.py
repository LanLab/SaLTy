from time import sleep as sl
import sys
import unittest
from unittest.mock import patch
from salty import salty
import os
import argparse
import shutil
import pandas as pd

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
        args = argparse.Namespace(input_folder = cwd+"/salty/tests/inputs/",
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
        sl(0.5)
        outsummary = self.outputfolder+"/summaryReport.txt"

        expected_outsummary = cwd+"/salty/tests/expected_outputs/summaryReport.txt"


        self.assertListEqual(
            list(open(outsummary)),
            list(open(expected_outsummary)))

    def tearDown(self):
        shutil.rmtree(self.outputfolder)

class TestCollectGenomes(unittest.TestCase):

    def test_collectGenomes(self):
        # Create a temporary folder and some test files
        test_folder = "test_folder"
        os.makedirs(test_folder, exist_ok=True)
        open(os.path.join(test_folder, "test.fasta"), 'w').close()
        open(os.path.join(test_folder, "test.fna"), 'w').close()
        open(os.path.join(test_folder, "test_1.fastq.gz"), 'w').close()
        open(os.path.join(test_folder, "test_2.fastq.gz"), 'w').close()

        # Create Namespace object simulating command-line arguments
        args = argparse.Namespace(input_folder=test_folder)

        # Call the function
        paths = salty.collectGenomes(args)

        # Assert that paths are collected correctly
        expected_paths = [
            ['assembly', os.path.join(test_folder, 'test.fasta')],
            ['assembly', os.path.join(test_folder, 'test.fna')],
            ['pairedEndReadForward', os.path.join(test_folder, 'test_1.fastq.gz')]
        ]
        self.assertEqual(paths, expected_paths)

        # Clean up: Remove temporary folder and files
        for file in os.listdir(test_folder):
            os.remove(os.path.join(test_folder, file))
        os.rmdir(test_folder)

class TestCheckInputReads(unittest.TestCase):

    def test_checkInputReads_with_valid_input(self):
        # Create a temporary folder and some test files
        test_folder = "test_folder"
        os.makedirs(test_folder, exist_ok=True)
        open(os.path.join(test_folder, "test_1.fastq.gz"), 'w').close()
        open(os.path.join(test_folder, "test_2.fastq.gz"), 'w').close()

        # Call the function
        result = salty.checkInputReads(test_folder)

        # Assert that the result is True when there are two 'fastq.gz' files
        self.assertTrue(result)

        # Clean up: Remove temporary folder and files
        for file in os.listdir(test_folder):
            os.remove(os.path.join(test_folder, file))
        os.rmdir(test_folder)

    def test_checkInputReads_with_invalid_input(self):
        # Create a temporary folder and some test files
        test_folder = "test_folder"
        os.makedirs(test_folder, exist_ok=True)
        open(os.path.join(test_folder, "test_1.fastq.gz"), 'w').close()

        # Call the function
        result = salty.checkInputReads(test_folder)

        # Assert that the result is False when there is only one 'fastq.gz' file
        self.assertFalse(result)

        # Clean up: Remove temporary folder and files
        for file in os.listdir(test_folder):
            os.remove(os.path.join(test_folder, file))
        os.rmdir(test_folder)

    def test_checkInputReads_with_empty_folder(self):
        # Create a temporary empty folder
        test_folder = "empty_folder"
        os.makedirs(test_folder, exist_ok=True)

        # Call the function
        result = salty.checkInputReads(test_folder)

        # Assert that the result is False when the folder is empty
        self.assertFalse(result)

        # Clean up: Remove temporary folder
        os.rmdir(test_folder)

# TODO troubleshoot test errors
# class TestGetLineageFromAllele(unittest.TestCase):
#     def setUp(self):
#         # Sample data for testing
#         self.alleles = pd.DataFrame({
#             'Allele': ['A', 'B', 'C'],
#             # Add more columns as needed for your test data
#         })
#         self.args = argparse.Namespace(lineages='lineages.csv')  # Sample args
#         self.outpath = 'output/'
# 
#     def test_single_lineage(self):
#         # Test when only one lineage is found
#         expected_lineage = 'LineageA'
#         expected_result = self.alleles.copy()
#         expected_result['Lineage'] = expected_lineage
# 
#         # Mock the pd.read_csv() function to return a DataFrame with a single lineage
#         with patch('pandas.read_csv', return_value=pd.DataFrame({'Lineage': [expected_lineage]})):
#             result = salty.getLineageFromAllele(self.alleles, self.args, self.outpath)
# 
#         self.assertTrue(expected_result.equals(result))
# 
#     def test_no_lineage(self):
#         # Test when no lineage is found
#         expected_result = self.alleles.copy()
#         expected_result['Lineage'] = 'No lineages association.'
# 
#         # Mock the pd.read_csv() function to return an empty DataFrame
#         with patch('pandas.read_csv', return_value=pd.DataFrame()):
#             result = salty.getLineageFromAllele(self.alleles, self.args, self.outpath)
# 
#         self.assertTrue(expected_result.equals(result))
# 
#     def test_multiple_lineages(self):
#         # Test when multiple lineages are found
#         expected_lineages = pd.DataFrame({
#             'Lineage': ['LineageA', 'LineageB']
#         })
#         expected_result = self.alleles.copy()
#         expected_result['Lineage'] = 'Mulitple Lineage Found'
# 
#         # Mock the pd.read_csv() function to return a DataFrame with multiple lineages
#         with patch('pandas.read_csv', return_value=expected_lineages):
#             result = salty.getLineageFromAllele(self.alleles, self.args, self.outpath)
# 
#         # Check if the output file is created
#         self.assertTrue(os.path.exists(self.outpath + '_multipleLineageAlleles.csv'))
#         self.assertTrue(expected_result.equals(result))

# TODO troubleshoot test errors
# class TestFiltLineageAlleles(unittest.TestCase):
#     def test_filtLineageAlleles(self):
#         # Sample data for testing
#         alleles = {'Lineage': '-', 'GeneA': 'A', 'GeneB': 'T'}
#         lineageAllelesDF = pd.DataFrame({
#             'GeneA': ['A', 'S', 'C'],
#             'GeneB': ['X', 'Y', 'Z'],
#             'Lineage': ['Lineage1', 'Lineage2', 'Lineage3']
#         })
# 
#         # Expected output
#         expected_result = pd.DataFrame({
#             'GeneA': ['A'],
#             'GeneB': ['T'],
#             'Lineage': ['Lineage1']
#         })
# 
#         result = salty.filtLineageAlleles(alleles, lineageAllelesDF)
#         resultc = list(sorted(result.columns))
#         result = result[resultc]
#         print(result)
#         print(expected_result)
#         # Check if the result matches the expected output
#         self.assertTrue(expected_result.equals(result))