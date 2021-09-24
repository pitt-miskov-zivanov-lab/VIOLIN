"""
in_out.py

Handles file input and output functions for VIOLIN tool
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd
import os.path
import numpy as np

from VIOLIN.formatting import add_regulator_names_id, evidence_score
from VIOLIN.network import node_edge_list


# Default Kind Score values
kind_dict = {"strong corroboration": 2,
             "weak corroboration1": 1,
             "weak corroboration2": 1,
             "weak corroboration3": 1,
             "hanging extension": 40,
             "full extension": 40,
             "internal extension": 40,
             "specification": 30,
             "dir contradiction": 10,
             "sign contradiction": 11,
             "att contradiction": 12,
             "flagged1": 20,
             "flagged2": 20,
             "flagged3": 20}


# Default Column Names for reading and model spreadsheets
reading_columns = ['Element Name', 'Element Type', 'Element ID',
                   'Positive Reg Name', 'Positive Reg Type', 'Positive Reg ID',
                   'Negative Reg Name', 'Negative Reg Type', 'Negative Reg ID',
                   'Connection Type', 'Mechanism', 'Paper ID', 'Evidence']
model_columns = ['Element Name', 'Element Type', 'Element IDs', 'Variable',
                 'Positive Regulators', 'Positive Regulators Connection Type',
                 'Negative Regulators', 'Negative Regulators Connection Type']

# Necessary column names for VIOLIN to work properly
required_model = ['Element Name', 'Element Type', 'Element IDs', 'Variable',
                  'Positive Regulators', 'Negative Regulators']

# Default Column names for calculating evidence score
evidence_score_def = ['Element Name','Element Type','Element ID','Positive Reg Name',
                       'Positive Reg Type','Positive Reg ID','Negative Reg Name',
                       'Negative Reg Type','Negative Reg ID','Connection Type']

def input_biorecipes(model, model_cols=model_columns):
    """
    This function imports a model file which is already in the BioRECIPES
    format, and converts all characters to lower case

    Parameters
    ----------
    model : str
        Directory and filename of the file containing the model spreadsheet
        in BioRECIPES format
        Accepted files: .txt, .csv, .tsv, .xlsx
    model_cols : list
        Column names of the model file. Default names found in model_columns

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

    # Check to make sure column names in the files are compatible with
    # the Accepted Column Names and Required Column Names
    if {(set(model_cols).issubset(set(model_df.columns))) and
        (set(required_model).issubset(set(model_cols)))}:
        # Remove extraaneous whitespace
        model_df = model_df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
        # Covert regulator variable name lists to common names
        # and database identifiers
        new_model = add_regulator_names_id(model_df)
        # Convert all model text to lower-case
        new_model = new_model.apply(lambda x: x.astype(str).str.lower())

        # Check that dataframes have been created correctly
        # new_model.to_csv(r'MODEL_dataframe.csv')
    else:
        raise ValueError("Either your file does not match the column names," +
                         " or you are missing necessary columns" + "\n" +
                         "The model file column names you input: " + str(model_cols) + "\n" +
                         "VIOLIN requires the following model columns: "+str(required_model))
    return new_model


def input_reading(reading, evidence_score_cols=evidence_score_def, atts=[]):
    """
    This function imports the reading file into the correct mode

    Parameters
    ----------
    reading : str
        Directory and filename of the machine reading spreadsheet output
        Accepted files: .txt, .csv, .tsv, .xlsx
    evidence_score_cols : list
        Column headings used to identify identical interactions in the machine reading output
    atts : list
        List of additional attributes which are available in LEE output
        Default is none

    Returns
    -------
    new_reading : pd.DataFrame]
        Formatted reading dataframe, including evidence count and list of PMCIDs
    """

    #Upload the model and reading files as dataframes based on the file extension
    reading_ext = os.path.splitext(reading)[1]

    if reading_ext == '.txt': reading_df = pd.read_csv(reading, sep='\t',index_col=None).fillna("nan")
    elif reading_ext == '.csv': reading_df = pd.read_csv(reading, sep=',',index_col=None).fillna("nan")
    elif reading_ext == '.xlsx': reading_df = pd.read_excel(reading,index_col=None).fillna("nan")
    elif reading_ext == '.tsv': reading_df = pd.read_csv(reading, sep='\t',index_col=None).fillna("nan")
    else: raise ValueError("The accepted file extensions are .txt, .csv, .xslx, and .tsv")

    #Begin relative column name retrieval 
    #Accepted target/regulated headers
    t_name_list = ["elementname","targetname","regulatedname"]
    t_type_list = ["elementtype","targettype","regulatedtype"]
    t_id_list = ["elementid","elementidentifier","targetid","targetidentifier","regulatedid","regulatredidentifier"]
    t_att_pre = ["element","target","regulated",""]

    #Accepted source/regulator headers (assuming separate positive and negative columns)
    s_name_list = ["positiveregname","posregname","positiveregulatorname","posregulatorname",
                    "positivesourcename","possourcename"]
    s_type_list = ["positiveregtype","posregtype","positiveregulatortype","positiveregtype",
                    "positivesourcetype","possourcetype"]
    s_id_list = ["positiveregid","posregid","positiveregulatorid","positiveregid",
                    "positivesourceid","possourceid","positiveregidentifier","posregidentifier",
                    "positiveregulatoridentifier","positiveregidentifier","positivesourceidentifier","possourceidentifier"]
    s_att_pre = ["posreg","posregulator","positivereg","positiveregulator","possource","positivesource"]

    #Get column header names as list
    col_names = list(reading_df.columns)

    #formate to bare minimum information
    bare_cols = [x.lower().replace(" ","").replace("_","").replace("-","") for x in col_names]

    #Check intersection of accepted column names and file column names
    if {len(set(t_name_list) & set(bare_cols)) == 1 & len(set(t_type_list) & set(bare_cols)) == 1 & 
        len(set(t_id_list) & set(bare_cols)) == 1 & len(set(s_name_list) & set(bare_cols)) == 1 & 
        len(set(s_type_list) & set(bare_cols)) == 1 & len(set(s_id_list) & set(bare_cols)) == 1}:
        #If minimum necessary columns are found, define variables for the column header
        target_name = col_names[bare_cols.index((set(t_name_list) & set(bare_cols)).pop())]
        target_type = col_names[bare_cols.index((set(t_type_list) & set(bare_cols)).pop())]
        target_id = col_names[bare_cols.index((set(t_id_list) & set(bare_cols)).pop())]
        pos_source_name = col_names[bare_cols.index((set(s_name_list) & set(bare_cols)).pop())]
        pos_source_type = col_names[bare_cols.index((set(s_type_list) & set(bare_cols)).pop())]
        pos_source_id = col_names[bare_cols.index((set(s_id_list) & set(bare_cols)).pop())]
        neg_source_name = pos_source_name.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
        neg_source_type = pos_source_type.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
        neg_source_id = pos_source_id.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
        #store column header names in a dictionary
        reading_cols = {"target_name" : target_name,"target_type" :target_type,"target_id" : target_id,
                        "pos_source_name" : pos_source_name, "pos_source_type" : pos_source_type, "pos_source_id" : pos_source_id,
                        "neg_source_name" : neg_source_name, "neg_source_type" : neg_source_type, "neg_source_id" : neg_source_id}
        #Now for the attributes:
        for x in atts:
            t_att_list = [pre + x.lower().replace(" ","") for pre in t_att_pre]
            s_att_list = [pre + x.lower().replace(" ","") for pre in s_att_pre]
            #made sure the attribute is in the reading columns
            if len(set(t_att_list) & set(bare_cols)) == 1 & len(set(s_att_list) & set(bare_cols)) == 1:
                #if it's found, add it to the reading columns
                #add the attribute for the target
                reading_cols['target_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(t_att_list) & set(bare_cols)).pop())]
                reading_cols['pos_source_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(s_att_list) & set(bare_cols)).pop())]
                reading_cols['neg_source_'+x.lower().replace(" ","_")] = reading_cols['pos_source_'+x.lower().replace(" ","_")].replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
            else:
                raise ValueError("Attribute \""+x+"\" was not found in your LEE input document."+"\n"+
                "Please check your file and try again")
        
    else:
        raise ValueError("Your LEE input is missing information."+"\n"+
        "VIOLIN requires the following information: Name, Type, and ID of target node and regulators")
    #End relative column name retrieval

    #Make sure evidence_cols match what is in the LEE input file
    if (set(evidence_score_cols).issubset(set(reading_df.columns))): 
        #Calculate the Evidence Score
        new_reading = evidence_score(reading_df,evidence_score_cols)
        ## Check that dataframes have been created correctly
        # new_reading.to_csv(r'READING_dataframe.csv')
    else: raise ValueError("The columns you chose for calculating the Evidence Score are not in youe LEE input file:"+str(evidence_score_cols))
    return new_reading, reading_cols


def output(reading_df, file_name, kind_values=kind_dict):
    """
    This function outputs the scored reading interactions.
    This writes output files, there are no return variables

    Parameters
    ----------
    reading_df : pd.DataFrame
        Dataframe of the scored reading dataframe
    file_name : str
        Directory and filename of the output suffix
    kind_values : dict
        Dictionary containing the numerical values for the Kind Score classifications
        Default values are found in kind_dict
    """

    #Output with all reading interactions, sorted by highest Total Score
    outputdf = reading_df.sort_values(by='Total Score', ascending=False)
    output_file = file_name+'_TotalOutput.csv'
    outputdf.to_csv(output_file)

    ## Corroborations ##
    corr = reading_df[(reading_df['Kind Score'] == kind_values['strong corroboration']) |
                      (reading_df['Kind Score'] == kind_values['weak corroboration1']) |
                      (reading_df['Kind Score'] == kind_values['weak corroboration2']) |
                      (reading_df['Kind Score'] == kind_values['weak corroboration3'])]
    corr = corr.sort_values(by='Total Score', ascending=False)
    corr_file = file_name+'_corroborations.csv'
    corr.to_csv(corr_file)

    ## Extensions ##
    ext = reading_df[(reading_df['Kind Score'] == kind_values['hanging extension']) | (reading_df['Kind Score'] == kind_values['full extension']) | (reading_df['Kind Score'] == kind_values['internal extension']) | (reading_df['Kind Score'] == kind_values['specification'])]
    ext = ext.sort_values(by='Total Score', ascending=False)
    ext_file = file_name+'_extensions.csv'
    ext.to_csv(ext_file)

    ## Contradictions ##
    cont = reading_df[(reading_df['Kind Score'] == kind_values['dir contradiction']) | (reading_df['Kind Score'] == kind_values['sign contradiction']) | (reading_df['Kind Score'] == kind_values['att contradiction'])]
    cont = cont.sort_values(by='Total Score', ascending=False)
    cont_file = file_name+'_contradictions.csv'
    cont.to_csv(cont_file)

    ## Special Cases ##
    que = reading_df[(reading_df['Kind Score'] == kind_values['flagged1']) | (reading_df['Kind Score'] == kind_values['flagged2']) | (reading_df['Kind Score'] == kind_values['flagged3'])]
    que = que.sort_values(by='Total Score', ascending=False)
    que_file = file_name+'_flagged.csv'
    que.to_csv(que_file)

    return
