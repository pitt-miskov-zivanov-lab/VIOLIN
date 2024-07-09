"""
in_out.py

Handles file input and output functions for VIOLIN tool
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd
import os.path
import numpy as np
import warnings
from violin.formatting import add_regulator_names_id, evidence_score, get_element, format_variable_names, wrap_list_to_str, get_listname
from violin.network import node_edge_list
import warnings
import re

# Default Kind Score values
kind_dict = {"strong corroboration": 2,
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



model_columns = ['Element Name', 'Element Type', 'Element IDs', 'Variable',
                 'Positive Regulator List', 'Positive Connection Type List',
                 'Negative Regulator List', 'Negative Connection Type List']

# Necessary column names for VIOLIN to work properly
required_model = ['Element Name', 'Element Type', 'Element IDs', 'Variable',
                  'Positive RegulatorList', 'Negative Regulator List']

# Default Column names for calculating evidence score
evidence_score_def = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol", "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                        "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol", "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                        "Sign", "Connection Type", "Mechanism", "Site",
                        "Cell Line", "Cell Type", "Tissue Type", "Organism"]

BioRECIPE_reading_col = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol", "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                        "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol", "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                        "Sign", "Connection Type", "Mechanism", "Site",
                        "Cell Line", "Cell Type", "Tissue Type", "Organism",
                        "Score", "Source", "Statements", "Paper IDs"]
def preprocessing_model(model, model_cols=model_columns):
    """
    This function check if your model is correct or necessary columns are missing or not
    Parameters
    ----------
    model : str
        model filename
    model_cols : list
        A list of model column names
    Returns
    -------
    new_model : pd.DataFrame
        Formatted model dataframe
    """
    # Upload the model and reading files as dataframes based on the file extension
    model_ext = os.path.splitext(model)[1]

    if model_ext == '.txt': model_df = pd.read_csv(model, sep='\t', index_col=None).fillna("nan")
    elif model_ext == '.csv': model_df = pd.read_csv(model, sep=',', index_col=None).fillna("nan")
    elif model_ext == '.xlsx': model_df = pd.read_excel(model, index_col=None).fillna("nan")
    elif model_ext == '.tsv': model_df = pd.read_csv(model, sep='\t',index_col=None).fillna("nan")
    else: raise ValueError("The accepted file extensions are .txt, .csv, .xslx, and .tsv")

    model_index = model_df.index
    model_df = format_variable_names(model_df)

    if {(set(model_cols).issubset(set(model_df.columns))) and
        (set(required_model).issubset(set(model_cols)))}:

    # FIXME: this block for getting list from regulation rule will be deprecated

    #     # Check if the model has the regulator information
    #     pos_regulator = all(x in ['', 'Nan', 'nan'] for x in model_df['Positive Regulator List'])
    #     neg_regulator = all(x in ['', 'Nan', 'nan'] for x in model_df['Negative Regulator List'])
    #     pos_regulation = all(x in ['', 'Nan', 'nan'] for x in model_df[f'Positive Regulation Rule'])
    #     neg_regulation = all(x in ['', 'Nan', 'nan'] for x in model_df[f'Negative Regulation Rule'])
    #
    #     if (pos_regulator and pos_regulation and neg_regulator and neg_regulation):
    #         raise ValueError(
    #             "The regulation rule and list columns are both empty, please fill at least one column out")
    #     elif (not pos_regulator) or (not neg_regulator):
    #         pass
    #     # Get regulator list from regulation rule if regulator list is empty
    #     else:
    #         for row in range(len(model_df)):
    #             index = model_index[row]
    #             for sign in ['Positive', 'Negative']:
    #                 if model_df.loc[index, f'{sign} Regulation Rule'] == '':
    #                     model_df.loc[index, f'{sign} Regulator List'] = ''
    #                 else:
    #                     model_df.loc[index, f'{sign} Regulator List'] = ','.join(
    #                         list(get_element(model_df.loc[index, f'{sign} Regulation Rule'], 0)))

        # Create a column for list-name
        model_df['Listname'] = [get_listname(idx, model_df) for idx in range(len(model_df))]
        # Normalize element type
        model_df['Element Type'] = model_df['Element Type'].str.replace(' ', '')
        # Covert regulator variable name lists to common names
        # and database identifiers
        new_model = add_regulator_names_id(model_df)
        # Remove extraaneous whitespace
        model_df = model_df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

        # Convert all model text to lower-case
        new_model = model_df.apply(lambda x: x.astype(str).str.lower())

    else:
        raise ValueError("Either your file does not match the column names," +
                             " or you are missing necessary columns" + "\n" +
                             "The model file column names you input: " + str(model_cols) + "\n" +
                             "VIOLIN requires the following model columns: " + str(required_model))
    return new_model


def preprocessing_reading(reading, evidence_score_cols=evidence_score_def, atts=[]):
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
    #Upload the model and reading files as dataframes based on the file extension
    reading_ext = os.path.splitext(reading)[1]

    if reading_ext == '.txt': reading_df = pd.read_csv(reading, sep='\t',index_col=None).fillna('nan')
    elif reading_ext == '.csv': reading_df = pd.read_csv(reading, sep=',',index_col=None).fillna('nan')
    elif reading_ext == '.xlsx': reading_df = pd.read_excel(reading,index_col=None).fillna('nan')
    elif reading_ext == '.tsv': reading_df = pd.read_csv(reading, sep='\t',index_col=None).fillna('nan')
    else: raise ValueError("The accepted file extensions are .txt, .csv, .xlsx, and .tsv")
    reading_df = reading_df.astype(str)
    for row in range(len(reading_df)):
        reading_df.loc[row, 'Regulator Type'] = ''.join(re.findall(r'[A-z]+', reading_df.loc[row, 'Regulator Type'].lower())) \
                                    if reading_df.loc[row, 'Regulator Type'] != 'nan' else reading_df.loc[row, 'Regulator Type']
        reading_df.loc[row, 'Regulated Type'] = ''.join(re.findall(r'[A-z]+', reading_df.loc[row, 'Regulated Type'].lower())) \
                                    if reading_df.loc[row, 'Regulated Type'] != 'nan' else reading_df.loc[row, 'Regulated Type']

        if reading_df.loc[row, 'Connection Type'].lower() in ['I', 'i', 'indirect', 'false']:
            reading_df.loc[row, 'Connection Type'] = 'i'
        elif reading_df.loc[row, 'Connection Type'].lower() in ['', 'nan', 'none']:
            reading_df.loc[row, 'Connection Type'] = 'i'
            warnings.warn(f'Connection type does not exist in row {row}, saving as indirect connection type.')
        else:
            reading_df.loc[row, 'Connection Type'] = 'd'
    #Make sure evidence_cols match what is in the LEE input file
    if (set(evidence_score_cols).issubset(set(reading_df.columns))):
        #Calculate the Evidence Score
        new_reading = evidence_score(reading_df,evidence_score_cols)
    else: raise ValueError("The columns you chose for calculating the Evidence Score are not in youe LEE input file:"+str(evidence_score_cols))
    return new_reading

def output(reading_df, file_name, kind_values=kind_dict):
    """
    This function outputs the scored reading interactions.
    This writes output files, there are no return variables

    Parameters
    ----------
    reading_df : pd.dataframe
        dataframe of the scored reading dataframe
    file_name : str
        Directory and filename of the output suffix
    kind_values : dict
        Dictionary containing the numerical values for the Kind Score classifications
        Default values are found in kind_dict
    """
    global BioRECIPE_reading_col

    #reading_df.reset_index(inplace=True)
    reading_df = reading_df.replace('nan', '')
    reading_df = wrap_list_to_str(reading_df, ['Score', 'Source', 'Statements', 'Paper IDs'])
    reading_df[BioRECIPE_reading_col] = reading_df[BioRECIPE_reading_col].astype(str)

    #Output with all reading interactions, sorted by highest Total Score
    outputdf = reading_df.sort_values(by='Total Score', ascending=False)
    outputdf.to_csv(f'{file_name}_outputDF.csv', index=False)
    output_file = file_name+'_scoreDF.csv'
    outputdf = outputdf[['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']]
    outputdf.to_csv(output_file, index=False)

    ## Corroborations ##
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

    ## Extensions ##
    ext = reading_df[(reading_df['Kind Score'] == kind_values['hanging extension']) |
                     (reading_df['Kind Score'] == kind_values['full extension']) |
                     (reading_df['Kind Score'] == kind_values['internal extension'])]
    ext = ext.sort_values(by='Total Score', ascending=False).reset_index()
    ext.to_csv(f'{file_name}_extensions.csv', index=False)
    output_file = file_name + '_extensions_score.csv'
    ext = ext[['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']]
    ext.to_csv(output_file, index=False)

    ## Contradictions ##
    cont = reading_df[(reading_df['Kind Score'] == kind_values['dir contradiction']) |
                      (reading_df['Kind Score'] == kind_values['sign contradiction']) |
                      (reading_df['Kind Score'] == kind_values['att contradiction'])]
    cont = cont.sort_values(by='Total Score', ascending=False).reset_index()
    cont.to_csv(f'{file_name}_contradictions.csv', index=False)
    output_file = file_name + '_contradictions_score.csv'
    cont = cont[['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']]
    cont.to_csv(output_file, index=False)

    ## Special Cases ##
    if ('flagged4' in kind_dict) and ('flagged5' in kind_dict):

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


# FIXME: This function will be deprecated
# def violin_to_biorecipe(VIOLIN_reading=None, VIOLIN_reading_df=None):
#     # FIXME: problematic right now, the classified output of VIOLIN is hard to be translated
#
#     """
#     Reverse function of BioRECIPE reading format to VIOLIN reading format
#     NOTE:
#     scoring columns will be added if VIOLIN_reading is the output
#     Returns
#     -------
#     BioRECIPE_reading: pd.DataFrame
#         Formatted reading dataframe in BioRECIPE format
#     """
#
#     if VIOLIN_reading:
#         VIOLIN_reading_df = pd.read_excel(VIOLIN_reading, index_col=None)
#
#     BioRECIPE_reading_df = pd.DataFrame(columns=BioRECIPE_reading_col)
#
#     # convert VIOLIN reading to BioRECIPE format
#     for i in range(len(VIOLIN_reading_df)):
#         if not VIOLIN_reading_df.loc[i, "Element Name"]:
#             break
#
#         BioRECIPE_reading_df.loc[i, "Regulated Name"] = VIOLIN_reading_df.loc[i, "Element Name"]
#         BioRECIPE_reading_df.loc[i, "Regulated Type"] = VIOLIN_reading_df.loc[i, "Element Type"]
#         BioRECIPE_reading_df.loc[i, "Regulated ID"] = VIOLIN_reading_df.loc[i, "Element ID"]
#         BioRECIPE_reading_df.loc[i, "Regulated Compartment"] = VIOLIN_reading_df.loc[i, "Location"]
#         BioRECIPE_reading_df.loc[i, "Regulated Compartment ID"] = VIOLIN_reading_df.loc[i, "Location ID"]
#
#         if VIOLIN_reading_df.loc[i, "Positive Reg Name"]:
#             BioRECIPE_reading_df.loc[i, "Regulator Name"] = VIOLIN_reading_df.loc[i, "Positive Reg Name"]
#             BioRECIPE_reading_df.loc[i, "Regulator Type"] = VIOLIN_reading_df.loc[i, "Positive Reg Type"]
#             BioRECIPE_reading_df.loc[i, "Regulator ID"] = VIOLIN_reading_df.loc[i, "Positive Reg ID"]
#             BioRECIPE_reading_df.loc[i, "Regulator Compartment"] = VIOLIN_reading_df.loc[i, "Positive Reg Location"]
#             BioRECIPE_reading_df.loc[i, "Regulator Compartment ID"] = VIOLIN_reading_df.loc[i, "Positive Reg Location ID"]
#             BioRECIPE_reading_df.loc[i, "Sign"] = "positive"
#         elif VIOLIN_reading_df.loc[i, "Negative Reg Name"]:
#             BioRECIPE_reading_df.loc[i, "Regulator Name"] = VIOLIN_reading_df.loc[i, "Negative Reg Name"]
#             BioRECIPE_reading_df.loc[i, "Regulator Type"] = VIOLIN_reading_df.loc[i, "Negative Reg Type"]
#             BioRECIPE_reading_df.loc[i, "Regulator ID"] = VIOLIN_reading_df.loc[i, "Negative Reg ID"]
#             BioRECIPE_reading_df.loc[i, "Regulator Compartment"] = VIOLIN_reading_df.loc[i, "Negative Reg Location"]
#             BioRECIPE_reading_df.loc[i, "Regulator Compartment ID"] = VIOLIN_reading_df.loc[i, "Negative Reg Location ID"]
#             BioRECIPE_reading_df.loc[i, "Sign"] = "negative"
#         else:
#             raise ValueError("Element {0} has no regulator".format(VIOLIN_reading_df.loc[i, "Element Name"]))
#
#         if VIOLIN_reading_df.loc[i, "Connection Type"].lower() == "d":
#             BioRECIPE_reading_df.loc[i, "Connection Type"] = "True"
#         elif VIOLIN_reading_df.loc[i, "Connection Type"].lower() == "i":
#             BioRECIPE_reading_df.loc[i, "Connection Type"] = "False"
#         else:
#             #print("Unspecified Connection Type: {0}".format(BioRECIPE_reading_df.loc[i, "Connection Type"]))
#             BioRECIPE_reading_df.loc[i, "Connection Type"] = "False"
#
#         BioRECIPE_reading_df.loc[i, "Mechanism"] = VIOLIN_reading_df.loc[i, "Mechanism"]
#         BioRECIPE_reading_df.loc[i, "Statements"] = VIOLIN_reading_df.loc[i, "Evidence"]
#         BioRECIPE_reading_df.loc[i, "Paper IDs"] = VIOLIN_reading_df.loc[i, "Paper ID"]
#
#         # TODO: specify this as VIOLIN's output with scoring columns
#         for score in ['Evidence Score', 'Match Score', 'Kind Score', 'Epistemic Value', 'Total Score']:
#             if VIOLIN_reading_df.loc[i, score]:
#                 BioRECIPE_reading_df.loc[i, score] = VIOLIN_reading_df.loc[i, score]
#
#     return BioRECIPE_reading_df