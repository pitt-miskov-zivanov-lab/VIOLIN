"""
use_Violin.py

Script used for calling full VIOLIN tool at the command line
Created November 2019 - Casey Hansen MeLoDy Lab
"""


import argparse
import os.path
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), os.pardir, '/src/violin')))

from violin.in_out import preprocessing_model, preprocessing_reading, output
from violin.scoring import score_reading
from violin.network import node_edge_list
from violin.visualize_violin import visualize

evidence_scoring_cols = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol", "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                        "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol", "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                        "Sign", "Connection Type", "Mechanism", "Site",
                        "Cell Line", "Cell Type", "Tissue Type", "Organism"]

attributes = ['Regulated Compartment ID', 'Regulator Compartment ID', 'Cell Line']

#Inputs: Model file, Reading File, Output Header, Classification, Filtering Option, Attributes
def use_violin(model_file, lee_file, out_file, approach = '1', score = 'extend', filt_opt = '100%', plot=True):
    """
    This function runs VIOLIN via a terminal command

    Parameters
    ----------
    model_file : str
        Directory and filename of the the machine reading spreadsheet output
        Accepted files: .txt, .csv, .tsv, .xlsx
    lee_file : str
        Directory and filename of the model file in BioRECIPE format
        Accepted files: .txt, .csv, .tsv, .xlsx
    out_file : str
         Directory and filename of the output suffix
    score : str
        Scoring scheme used for classification
        Options are: 'extend', 'extend subcategories', 'corroborate', 'corroborate subcategories'
    filt_opt : str
        How much VIOLIN output should be visualized. Can be filtered
        by top % of total score, evidence score (Se) threshold, or
        total score (St) threshold
        Accepted options are 'X%','Se>Y', or 'St>Z',
        where X, Y, and Z, are values
        Default is '100%' (Total Output)
    """
    # Defining the scoring scheme
    if score == 'extend':
        kind_dict = {"strong corroboration": 2,
                     "empty attribute": 1,
                     "indirect interaction": 3,
                     "path corroboration": 5,
                     "specification": 7,
                     "hanging extension": 40,
                     "full extension": 39,
                     "internal extension": 38,
                     "dir contradiction": 11,
                     "sign contradiction": 10,
                     "att contradiction": 9,
                     "dir mismatch": 20,
                     "path mismatch": 19,
                     "self-regulation": 18}
        match_dict = {"source present" : 1,
                        "target present" : 100,
                        "both present" : 10,
                        "neither present" : 0.1}
        if approach in ['1', '2']:
            pass
        else:
            kind_dict["flagged4"] = 20
            kind_dict["flagged5"] = 20


    elif score == 'corroborate':
        kind_dict = {"strong corroboration": 2,
                     "empty attribute": 1,
                     "indirect interaction": 3,
                     "path corroboration": 5,
                     "specification": 7,
                     "hanging extension": 40,
                     "full extension": 39,
                     "internal extension": 38,
                     "dir contradiction": 11,
                     "sign contradiction": 10,
                     "att contradiction": 9,
                     "dir mismatch": 20,
                     "path mismatch": 19,
                     "self-regulation": 18}
        match_dict = {"source present" : 1,
                        "target present" : 1,
                        "both present" : 100,
                        "neither present" : 0.1}

        if approach in ['1', '2']:
            pass
        else:
            kind_dict["flagged4"] = 1
            kind_dict["flagged5"] = 1

    elif score == 'extend subcategories':
        kind_dict = {"strong corroboration": 2,
                     "empty attribute": 1,
                     "indirect interaction": 3,
                     "path corroboration": 5,
                     "specification": 7,
                     "hanging extension": 40,
                     "full extension": 39,
                     "internal extension": 38,
                     "dir contradiction": 11,
                     "sign contradiction": 10,
                     "att contradiction": 9,
                     "dir mismatch": 20,
                     "path mismatch": 19,
                     "self-regulation": 18}
        match_dict = {"source present" : 1,
                        "target present" : 100,
                        "both present" : 10,
                        "neither present" : 0.1}

        if approach in ['1', '2']:
            pass
        else:
            kind_dict["flagged4"] = 23
            kind_dict["flagged5"] = 24

    elif score == 'corroborate subcategories':
        kind_dict = {"strong corroboration" : 40,
                        "weak corroboration1" : 30,
                        "weak corroboration2" : 31,
                        "weak corroboration3" : 32,
                        "hanging extension" : 2,
                        "full extension" : 4,
                        "internal extension" : 10,
                        "specification" : 11,
                        "dir contradiction" : 20,
                        "sign contradiction" : 21,
                        "att contradiction" : 22,
                        "flagged1" : 1,
                        "flagged2" : 3,
                        "flagged3" : 5}
        match_dict = {"source present" : 1,
                        "target present" : 1,
                        "both present" : 100,
                        "neither present" : 0.1}

        if approach in ['1', '2']:
            pass
        else:
            kind_dict["flagged4"] = 7
            kind_dict["flagged5"] = 9

    else:
        raise ValueError('Unaccepted scoring option'+'\n'+
                         'options are: \'extend\', \'extend subcategories\', \'corroborate\', \'corroborate subcategories\'')

    # Import model and LEE set, using default input parameters
    model_df = preprocessing_model(model_file)
    reading_df = preprocessing_reading(reading=lee_file,evidence_score_cols=evidence_scoring_cols, atts = attributes)
    graph = node_edge_list(model_df)

    #Scoring and Output
    scored = score_reading(reading_df,
                           model_df,
                           graph,
                           kind_values = kind_dict,
                           match_values = match_dict,
                           attributes=attributes,
                           classify_scheme = approach)
    output(scored,out_file,kind_values=kind_dict)

    #Visualization
    if plot:
        visualize(match_dict, kind_dict, out_file+'_outputDF.csv', filter_opt=filt_opt)
    else:
        pass

    return

def main():
    parser = argparse.ArgumentParser(description='Verifying Interactions Of Likely Importance to the Network')
    parser.add_argument('model', type=str,
                        help='file containing model interactions - must be extension .txt, .csv, .tsv, or .xlsx')
    parser.add_argument('reading', type=str,
                        help='file containing model interactions - must be extension .txt, .csv, .tsv, or .xlsx')
    parser.add_argument('output', type=str,
                        help='directory and suffix for output \n'
                        'example: /Users/casey/Desktop/PPC')
    parser.add_argument('score', type=str,
                        help='scoring value goal: extention, extension subcategories, corroboration, corroboration categories')
    parser.add_argument('filter', type=str,
                        help='(optional) filtering value for visualization output')

    parser.add_argument('approach', type=str, choices=['1', '2', '3'],
                        help='(optional) classify schemes, default is 1')
    args = parser.parse_args()

    if (os.path.splitext(args.model)[1] in ['.txt','.csv','.tsv','.xlsx'] and os.path.splitext(args.reading)[1] in ['.txt','.csv','.tsv','.xlsx'] and type(args.output)==str):
        if args.filter == None:
            if args.approach == None:
                use_violin(args.model,args.reading,args.output,args.score)
            else:
                use_violin(args.model,args.reading,args.output,args.approach,args.score)
        else:
            if args.approach == None:
                use_violin(args.model,args.reading,args.output,args.score,args.filter)
            else:
                use_violin(args.model,args.reading,args.output,args.approach,args.score,args.filter)

    else:
        raise ValueError('Unrecognized input format')
        # print('Unrecognized input format',os.path.splitext(args.reading)[1])

if __name__ == '__main__':
    main()
