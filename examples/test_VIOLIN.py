## This is script for testing the functionality of VIOLIN
import unittest
import pandas as pd
import os
import sys

# sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), os.pardir, 'src/violin/')))
# from in_out import preprocessing_model, preprocessing_reading, output
# from scoring import score_reading
# from network import node_edge_list
# from visualize_violin import visualize

from use_violin_script import use_violin

kind_dict = {"strong corroboration" : 2,
                "empty attribute" : 1,
                "indirect interaction" : 3,
                "path corroboration" : 5,
                "specification" : 7,
                "hanging extension" : 40,
                "full extension" : 39,
                "internal extension" : 38,
                "dir contradiction" : 11,
                "sign contradiction" : 10,
                "att contradiction" : 9,
                "dir mismatch" : 20,
                "path mismatch" : 19,
                "self-regulation" : 18}

match_dict = {"source present" : 1,
                "target present" : 100,
                "both present" : 10,
                "neither present" : 0.1}

evidence_scoring_cols = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol", "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                        "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol", "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                        "Sign", "Connection Type", "Mechanism", "Site",
                        "Cell Line", "Cell Type", "Tissue Type", "Organism"]

attributes = ['Regulated Compartment ID', 'Regulator Compartment ID', 'Cell Line']

model_file = 'input/models/SkMel133_biorecipe.xlsx'

class TestVIOLIN(unittest.TestCase):

    # Test Corroborations
    def test_strong_corroboration(self):
        reading_file = 'test/input/strong_corroboration.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_corroborations.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['strong corroboration'])

    def test_empty_attribute(self):
        reading_file = 'test/input/empty_attribute.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_corroborations.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['empty attribute'])

    def test_indirect_interaction(self):
        reading_file = 'test/input/indirect_interaction.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_corroborations.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['indirect interaction'])

    def test_path_corroboration(self):
        reading_file = 'test/input/path_corroboration.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_corroborations.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['path corroboration'])

    def test_specification(self):
        reading_file = 'test/input/specification.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_corroborations.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['specification'])

    # Test contradictions
    def test_dir_contradiction(self):
        reading_file = 'test/input/dir_contradiction.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_contradictions.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['dir contradiction'])

    def test_att_contradiction(self):
        reading_file = 'test/input/att_contradiction.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_contradictions.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['att contradiction'])

    def test_sign_contradiction(self):
        reading_file = 'test/input/sign_contradiction.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_contradictions.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['sign contradiction'])

    # Test flagged interactions
    # FIXME: dir mismatch are not able to test since dir contradiction is prior to it.
    # def test_dir_mismatch(self):
    #     reading_file = 'test/input/dir_mismatch.xlsx'
    #     output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
    #     use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
    #     df = pd.read_csv(output_dir + '_flagged.csv', index_col=None)
    #     self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['dir mismatch'])

    def test_path_mismatch(self):
        reading_file = 'test/input/path_mismatch.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_flagged.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['path mismatch'])

    def test_self_regulation(self):
        reading_file = 'test/input/self_regulation.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_flagged.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['self-regulation'])

    # Test extensions
    def test_full_extension(self):
        reading_file = 'test/input/full_extension.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_extensions.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['full extension'])

    def test_hanging_extension(self):
        reading_file = 'test/input/hanging_extension.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_extensions.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['hanging extension'])

    def test_internal_extension(self):
        reading_file = 'test/input/internal_extension.xlsx'
        output_dir = 'test/output/' + reading_file.split('/')[-1].split('.')[0]
        use_violin(model_file, reading_file, output_dir, score='extend subcategories',plot=False)
        df = pd.read_csv(output_dir + '_extensions.csv', index_col=None)
        self.assertEqual(df.loc[0, 'Kind Score'], kind_dict['internal extension'])


if __name__ == '__main__':
    unittest.main()