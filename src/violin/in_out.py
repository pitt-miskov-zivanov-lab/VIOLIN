"""
in_out.py

Handles file input and output functions for VIOLIN tool
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import os.path
import warnings

import pandas as pd
import re

from violin.formatting import get_hgnc_symbol
import json

from violin.formatting import (
    evidence_score, get_type,
    format_variable_names, wrap_list_to_str, get_listname)

# Default Kind Score values
KIND_DICT = {"strong corroboration": 2,
             "empty attribute": 1,
             "indirect interaction": 1,
             "path corroboration": 1,
             "hanging extension": 40,
             "full extension": 40,
             "internal extension": 40,
             "specification": 30,
             "dir contradiction": 10,
             "sign contradiction": 11,
             "att contradiction": 12,
             "dir mismatch": 20,
             "path mismatch": 20,
             "self-regulation": 20}

MODEL_COLUMNS = ['Element Name', 'Element Type', 'Element IDs', 'Variable',
                 'Positive Regulator List', 'Positive Connection Type List',
                 'Negative Regulator List', 'Negative Connection Type List']

# Necessary column names for VIOLIN to work properly
REQUIRED_MODEL = ['Element Name', 'Element Type', 'Element IDs', 'Variable',
                  'Positive RegulatorList', 'Negative Regulator List']

# Default Column names for calculating evidence score
EVIDENCE_SCORE_DEF = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol",
                      "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                      "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol",
                      "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                      "Sign", "Connection Type", "Mechanism", "Site",
                      "Cell Line", "Cell Type", "Tissue Type", "Organism"]

BioRECIPE_READING_COL = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol",
                         "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                         "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol",
                         "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                         "Sign", "Connection Type", "Mechanism", "Site",
                         "Cell Line", "Cell Type", "Tissue Type", "Organism",
                         "Score", "Source", "Statements", "Paper IDs"]


def preprocessing_model(model):
    """
    This function check if your model is correct or necessary columns are missing or not

    Parameters
    ----------
    model : str
        model filename
    Returns
    -------
    new_model : pd.DataFrame
        Formatted model dataframe
    """
    # Upload the model and reading files as dataframes based on the file extension
    global MODEL_COLUMNS
    model_cols = MODEL_COLUMNS
    model_ext = os.path.splitext(model)[1]

    if model_ext == '.txt':
        model_df = pd.read_csv(model, sep='\t', index_col=None).fillna("nan")
    elif model_ext == '.csv':
        model_df = pd.read_csv(model, sep=',', index_col=None).fillna("nan")
    elif model_ext == '.xlsx':
        model_df = pd.read_excel(model, index_col=None).fillna("nan")
    elif model_ext == '.tsv':
        model_df = pd.read_csv(model, sep='\t', index_col=None).fillna("nan")
    else:
        raise ValueError("The accepted file extensions are .txt, .csv, .xslx, and .tsv")

    model_df = format_variable_names(model_df)

    if {(set(model_cols).issubset(set(model_df.columns))) and
        (set(REQUIRED_MODEL).issubset(set(model_cols)))}:

        # Create a column for list-name
        # model_df['Listname'] = [get_listname(idx, model_df) for idx in range(len(model_df))]
        model_df['Listname'] = model_df['Element Name'].to_list()
        # Normalize element type
        model_df['Element Type'] = model_df['Element Type'].str.replace(' ', '')
        # Covert regulator variable name lists to common names
        # and database identifiers
        # new_model = add_regulator_names_id(model_df)
        # Remove extraaneous whitespace
        model_df = model_df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

        # Convert all model text to lower-case
        new_model = model_df.apply(lambda x: x.astype(str).str.lower())
        new_model['Element Type'] = new_model['Element Type'].apply(lambda x: get_type(x))  # normalize type

    else:
        raise ValueError("Either your file does not match the column names," +
                         " or you are missing necessary columns" + "\n" +
                         "The model file column names you input: " + str(model_cols) + "\n" +
                         "VIOLIN requires the following model columns: " + str(REQUIRED_MODEL))
    return new_model


def preprocessing_reading(reading, evidence_score_cols=None, atts=None):
    """
    This function import the reading file and check if the reading format is correct

    Parameters
    ----------
    reading : str
        Directory and filename of the machine reading spreadsheet output, in BioRECIPE format
        Accepted file: .txt, .csv, .tsv, .xlsx
    evidence_score_cols : list
        Column headings used to identify identical interactions in the machine reading output
    atts : list
        a List of additional attributes which are available in LEE output
        Default is none

    Returns
    -------
    new_reading : pd.dataframe
        formatted reading dataframe, including evidence count and list of PMCIDs
    """
    # Upload the model and reading files as dataframes based on the file extension
    # initialize default values in function
    if evidence_score_cols is None:
        evidence_score_cols = EVIDENCE_SCORE_DEF

    if atts is None:
        atts = []
    reading_ext = os.path.splitext(reading)[1]

    read_functions = {
        '.txt': pd.read_csv,
        '.csv': pd.read_csv,
        '.xlsx': pd.read_excel,
        '.tsv': pd.read_csv
    }

    if reading_ext not in read_functions:
        raise ValueError("The accepted file extensions are .txt, .csv, .xlsx, and .tsv")

    read_func = read_functions[reading_ext]
    kwargs = {'sep': '\t'} if reading_ext in ['.txt', '.tsv'] else {}
    reading_df = read_func(reading, index_col=None, **kwargs).fillna('nan')
    reading_df = reading_df.astype(str)
    reading_df = reading_df.applymap(lambda x: x.lower() if isinstance(x, str) else x)

    with open('hgnc_symbol_dict.json', 'r') as f:
        hgnc_dict = json.load(f)

    for row in range(len(reading_df)):

        if bool(re.fullmatch(r'[0-9\.]+', reading_df.loc[row, 'Regulator HGNC Symbol'])):
            try:
                reading_df.loc[row, 'Regulator HGNC Symbol'] = get_hgnc_symbol(reading_df.loc[row, 'Regulator HGNC Symbol'], hgnc_dict=hgnc_dict)
            except Exception as e:
                print(e)
        else:
            pass


        if bool(re.fullmatch(r'[0-9\.]+', reading_df.loc[row, 'Regulator HGNC Symbol'])):
            try:
                reading_df.loc[row, 'Regulated HGNC Symbol'] = get_hgnc_symbol(reading_df.loc[row, 'Regulated HGNC Symbol'], hgnc_dict=hgnc_dict)
            except Exception as e:
                print(e)
                pass                 
        else:
            pass

        if reading_df.loc[row, 'Connection Type'].lower() in ['I', 'i', 'indirect', 'false']:
            reading_df.loc[row, 'Connection Type'] = 'i'
        elif reading_df.loc[row, 'Connection Type'].lower() in ['', 'nan', 'none']:
            reading_df.loc[row, 'Connection Type'] = 'i'
            warnings.warn(f'Connection type does not exist in row {row}, saving as indirect connection type.')
        else:
            reading_df.loc[row, 'Connection Type'] = 'd'
    
    with open('hgnc_symbol_dict.json', 'w') as f:
        json.dump(hgnc_dict, f, indent=4)

    reading_df['Regulator Type'] = reading_df['Regulator Type'].apply(lambda x: get_type(x.lower()))  # normalize type
    reading_df['Regulated Type'] = reading_df['Regulated Type'].apply(lambda x: get_type(x.lower()))  # normalize type

    # Make sure evidence_cols match what is in the LEE input file
    if set(evidence_score_cols).issubset(set(reading_df.columns)):
        # Calculate the Evidence Score
        new_reading = evidence_score(reading_df, evidence_score_cols)
    else:
        raise ValueError(
            "The columns you chose for calculating the Evidence Score are not in youe LEE input file:" + str(
                evidence_score_cols))
    return new_reading


def output(reading_df, file_name, kind_values=None):
    """
    This function outputs the scored reading interactions.
    This writes output files, there are no return variables

    Parameters
    ----------
    reading_df : pd.dataframe
        the scored reading dataframe
    file_name : str
        Directory and filename of the output suffix
    kind_values : dict
        Dictionary containing the numerical values for the Kind Score classifications
        Default values are found in KIND_DICT
    """
    global BioRECIPE_READING_COL
    if kind_values is None:
        kind_values = KIND_DICT
    # reading_df.reset_index(inplace=True)
    reading_df = reading_df.replace('nan', '')
    reading_df = wrap_list_to_str(reading_df, ['Score', 'Source', 'Statements', 'Paper IDs'])
    reading_df[BioRECIPE_READING_COL] = reading_df[BioRECIPE_READING_COL].astype(str)

    # Output with all reading interactions, sorted by highest Total Score
    outputdf = reading_df.sort_values(by='Total Score', ascending=False)
    outputdf.to_csv(f'{file_name}_outputDF.csv', index=False)
    output_file = file_name + '_scoreDF.csv'
    outputdf = outputdf[['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']]
    outputdf.to_csv(output_file, index=False)

    # Corroborations #
    corr = reading_df[(reading_df['Kind Score'] == kind_values['strong corroboration']) |
                      (reading_df['Kind Score'] == kind_values['empty attribute']) |
                      (reading_df['Kind Score'] == kind_values['indirect interaction']) |
                      (reading_df['Kind Score'] == kind_values['path corroboration']) |
                      (reading_df['Kind Score'] == kind_values['specification'])]
    corr = corr.sort_values(by='Total Score', ascending=False).reset_index()
    corr.to_csv(f'{file_name}_corroborations.csv', index=False)
    output_file = file_name + '_corroborations_score.csv'
    corr = corr[['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']]
    corr.to_csv(output_file, index=False)

    # Extensions #
    ext = reading_df[(reading_df['Kind Score'] == kind_values['hanging extension']) |
                     (reading_df['Kind Score'] == kind_values['full extension']) |
                     (reading_df['Kind Score'] == kind_values['internal extension'])]
    ext = ext.sort_values(by='Total Score', ascending=False).reset_index()
    ext.to_csv(f'{file_name}_extensions.csv', index=False)
    output_file = file_name + '_extensions_score.csv'
    ext = ext[['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']]
    ext.to_csv(output_file, index=False)

    # Contradictions #
    cont = reading_df[(reading_df['Kind Score'] == kind_values['dir contradiction']) |
                      (reading_df['Kind Score'] == kind_values['sign contradiction']) |
                      (reading_df['Kind Score'] == kind_values['att contradiction'])]
    cont = cont.sort_values(by='Total Score', ascending=False).reset_index()
    cont.to_csv(f'{file_name}_contradictions.csv', index=False)
    output_file = file_name + '_contradictions_score.csv'
    cont = cont[['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']]
    cont.to_csv(output_file, index=False)

    # Special Cases #
    if ('flagged4' in KIND_DICT) and ('flagged5' in KIND_DICT):

        que = reading_df[(reading_df['Kind Score'] == kind_values['dir mismatch']) |
                         (reading_df['Kind Score'] == kind_values['path mismatch']) |
                         (reading_df['Kind Score'] == kind_values['self-regulation']) |
                         (reading_df['Kind Score'] == kind_values['flagged4']) |
                         (reading_df['Kind Score'] == kind_values['flagged5'])]
    else:
        que = reading_df[(reading_df['Kind Score'] == kind_values['dir mismatch']) |
                         (reading_df['Kind Score'] == kind_values['path mismatch']) |
                         (reading_df['Kind Score'] == kind_values['self-regulation'])]

    que = que.sort_values(by='Total Score', ascending=False).reset_index()
    que.to_csv(f'{file_name}_flagged.csv', index=False)
    output_file = file_name + '_flagged_score.csv'
    cont = cont[['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']]
    cont.to_csv(output_file, index=False)
    return
