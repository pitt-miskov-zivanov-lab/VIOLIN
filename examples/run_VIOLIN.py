"""
use_Violin.py

Script used for calling full VIOLIN tool at the command line
"""


import argparse
import os.path
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), os.pardir, '/src/violin')))

from violin.in_out import preprocessing_model, preprocessing_reading, output
from violin.scoring import score_reading
from violin.network import node_edge_list
from violin.visualize_violin import ViolinPlot

evidence_scoring_cols = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol", "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                        "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol", "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                        "Sign", "Connection Type", "Mechanism", "Site",
                        "Cell Line", "Cell Type", "Tissue Type", "Organism"]

attributes = ['Regulated Compartment ID', 'Regulator Compartment ID', 'Cell Line']

#Inputs: Model file, Reading File, Output Header, Classification, Filtering Option, Attributes
def use_violin(model_file, lee_file, out_file, approach = '1', attributes=attributes, score = 'extend', filt_opt = '100%', plot=True):
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
        Options are: \'corroboration\', \'extension\', \'contradiction\', \'flagged\', \'pie_plots\'
    filt_opt : str
        How much VIOLIN output should be visualized. Can be filtered
        by top % of total score, evidence score (Se) threshold, or
        total score (St) threshold
        Accepted options are 'X%','Se>Y', or 'St>Z',
        where X, Y, and Z, are values
        Default is '100%' (Total Output)
    """
    # Import model and LEE set, using default input parameters
    model_df = preprocessing_model(model_file)
    reading_df = preprocessing_reading(reading=lee_file,evidence_score_cols=evidence_scoring_cols, atts = attributes)
    graph = node_edge_list(model_df)

    #Scoring and Output
    scored = score_reading(reading_df,
                           model_df,
                           graph,
                           attributes=attributes,
                           classify_scheme = approach)
    output(scored, out_file, classify_scheme=approach)

    #Visualization
    print(plot)
    if plot:
        violin_plot = ViolinPlot(file_name=out_file+'_outputDF.csv', filter_opt=filt_opt, classify_scheme=approach)
        if score not in ['corroboration', 'extension', 'contradiction', 'flagged', 'pie_plots']:
            raise ValueError('Unrecognized scoring scheme')
        if score == 'pie_plots': 
            violin_plot.get_pie_plots(out_file=out_file, save=True, show=plot)
        else:
            violin_plot.get_category_summary(category=score, save_name=out_file+'_'+score, save=plot)
    else:
        pass

    return

def main():
    parser = argparse.ArgumentParser(description='Verifying Interactions Of Likely Importance to the Network')
    parser.add_argument('--model', type=str,
                        help='file containing model interactions - must be extension .txt, .csv, .tsv, or .xlsx')
    parser.add_argument('--reading', type=str,
                        help='file containing model interactions - must be extension .txt, .csv, .tsv, or .xlsx')
    parser.add_argument('--output', type=str,
                        help='directory and suffix for output \n'
                        'example: /Users/username/Desktop/PPC')
    parser.add_argument('--summary', type=str, default='corroboration',
                        help='(optional) the scoring scheme used for visualization.'\
                        'Options are: \'corroboration\', \'extension\', \'contradiction\', ' \
                        '\'flagged\', \'pie_plots\'')
    parser.add_argument('--filter', type=str, default='100%',
                        help='(optional) filtering value for visualization output')
    parser.add_argument('--attributes', type=str, nargs='+', default=attributes,
                        help='(optional) attributes to use for scoring and classification, default is ' \
                        'Regulated Compartment ID, Regulator Compartment ID, and Cell Line'\
                        'NOTE: Compartment and Compartment ID will be used for comparsion, as long as one of '\
                        'them is included in the list of attributes')
    parser.add_argument('--scheme', type=str, choices=['1', '2', '3'], default='1', 
                        help='(optional) classify schemes, default is 1')
    parser.add_argument('--show_plot', action='store_true',
                        help='(optional) whether to show the plot after the results are generated. Default is False')
    args = parser.parse_args()

    if (os.path.splitext(args.model)[1] in ['.txt','.csv','.tsv','.xlsx'] and os.path.splitext(args.reading)[1] in ['.txt','.csv','.tsv','.xlsx'] and type(args.output)==str):
        use_violin(model_file=args.model, 
                   lee_file=args.reading,
                   out_file=args.output,
                   approach=args.scheme,
                   attributes=args.attributes,
                   score=args.summary,
                   filt_opt=args.filter,
                   plot=args.show_plot)

    else:
        raise ValueError('Unrecognized input format')
        # print('Unrecognized input format',os.path.splitext(args.reading)[1])

if __name__ == '__main__':
    main()
