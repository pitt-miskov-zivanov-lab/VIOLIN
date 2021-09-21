"""
use_Violin.py

Script used for calling full VIOLIN tool at the command line
Created November 2019 - Casey Hansen MeLoDy Lab
"""


import argparse
import os.path

from VIOLIN.in_out import input_biorecipes,input_reading,output
from VIOLIN.scoring import score_reading
from VIOLIN.network import node_edge_list
from VIOLIN.visualize_violin import visualize

#Inputs: Model file, Reading File, Output Header, Classification, Filtering Option, Attributes
def use_violin(model_file, lee_file, out_file, score = 'extend', filt_opt = '100%'):
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
        kind_dict = {"strong corroboration" : 2, 
                    "weak corroboration1" : 1,
                    "weak corroboration2" : 1,
                    "weak corroboration3" : 1,
                    "hanging extension" : 40, 
                    "full extension" : 40, 
                    "internal extension" : 40, 
                    "specification" : 30, 
                    "dir contradiction" : 10,
                    "sign contradiction" : 10,
                    "att contradiction" : 10,
                    "flagged1" : 20,
                    "flagged2" : 20,
                    "flagged3" : 20}
        match_dict = {"source present" : 1, 
                        "target present" : 100, 
                        "both present" : 10, 
                        "neither present" : 0.1}
    elif score == 'corroborate':
        kind_dict = {"strong corroboration" : 40, 
                        "weak corroboration1" : 30,
                        "weak corroboration2" : 30,
                        "weak corroboration3" : 30,
                        "hanging extension" : 1, 
                        "full extension" : 1, 
                        "internal extension" : 10, 
                        "specification" : 10, 
                        "dir contradiction" : 20,
                        "sign contradiction" : 20,
                        "att contradiction" : 20,
                        "flagged1" : 1,
                        "flagged2" : 1,
                        "flagged3" : 1}
        match_dict = {"source present" : 1, 
                        "target present" : 1, 
                        "both present" : 100, 
                        "neither present" : 0.1}
    elif score == 'extend subcategories':
        kind_dict = {"strong corroboration" : 2, 
                    "weak corroboration1" : 1,
                    "weak corroboration2" : 3,
                    "weak corroboration3" : 5,
                    "hanging extension" : 40, 
                    "full extension" : 41, 
                    "internal extension" : 42, 
                    "specification" : 30, 
                    "dir contradiction" : 10,
                    "sign contradiction" : 11,
                    "att contradiction" : 12,
                    "flagged1" : 20,
                    "flagged2" : 21,
                    "flagged3" : 22}
        match_dict = {"source present" : 1, 
                        "target present" : 100, 
                        "both present" : 10, 
                        "neither present" : 0.1}
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
    else: 
        raise ValueError('Unaccepted scoring option'+'\n'+
                         'options are: \'extend\', \'extend subcategories\', \'corroborate\', \'corroborate subcategories\'')

    # Import model and LEE set, using default input parameters
    model_df = input_biorecipes(model_file)
    reading_df, reading_cols = input_reading(lee_file)
    graph = node_edge_list(model_df)

    #Scoring and Output
    scored = score_reading(reading_df,model_df,graph,reading_cols,kind_values = kind_dict,match_values = match_dict)
    output(scored,out_file,kind_values=kind_dict)
    
    #Visualization
    visualize(match_dict, kind_dict, out_file+'_TotalOutput.csv', filter_opt=filt_opt)

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
    

    args = parser.parse_args()

    if (os.path.splitext(args.model)[1] in ['.txt','.csv','.tsv','.xlsx'] and os.path.splitext(args.reading)[1] in ['.txt','.csv','.tsv','.xlsx'] and type(args.output)==str):
        if args.filter == None:
            use_violin(args.model,args.reading,args.output,args.score)
        else:
            use_violin(args.model,args.reading,args.output,args.score,args.filter)
    else: 
        raise ValueError('Unrecognized input format')
        # print('Unrecognized input format',os.path.splitext(args.reading)[1])

if __name__ == '__main__':
    main()
