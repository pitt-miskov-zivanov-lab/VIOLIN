import pandas as pd
import glob
import re
import sys
import os
from ast import literal_eval

sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), os.pardir, 'src/violin/')))
from numeric import find_element
from in_out import preprocessing_model

FILES_TEST = ["RB2"]
FILES = ["RA1",
         "RA2",
         "RA3",
         "RA4",
         "RB1",
         "RB2",
         "RB3",
         "RB_star_1",
         "RB_star_2"]

# FILES = ["RA1",
#          "RA2",
#          "RA3",
#          "RA4",
#          "RB1",
#          "RB2",
#          "RB3",
#          "RB_star_1",
#          "RB_star_2",
#          "RA2_0_1",
#          "RA2_0_1_1",
#          "RB2_1",
#          "RB2_0_1"]

def merge_duplicates(reading_df, col_names):
    """
    This function merges duplicate interactions and calculates evidence score of each LEE

    Parameters
    ----------
    reading_df : pd.DataFrame
        The dataframe of the machine reading output
    col_names: list
        Specifically the column headings used to determine if interactions are identical

    Returns
    -------
    counted_reading : pd.DataFrame
        A new dataframe with the evidence count and PMCID list for each interaction
    """

    #Convert reading to lower case, to prevent issues with case difference
    reading = reading_df.apply(lambda x: x.astype(str).str.lower())
    #The columns that aren't used to determine duplicates (such as Paper ID or Evidence Text)
    remainder = [x for x in reading_df.columns if x not in col_names]

    #As VIOLIN Identifies duplicates, it merges attributes from the remainder list into a single cell
    #This is how we count the number of times an LEE appears, and keep track of paper IDs and evidence text
    counted_reading = reading.groupby(col_names)[remainder[0]].apply(list).reset_index(name=remainder[0])
    for x in range(1,len(remainder)):
        sub = reading.groupby(col_names)[remainder[x]].apply(list).reset_index(name=remainder[x])
        counted_reading[remainder[x]] = sub[remainder[x]]

    return counted_reading

if __name__ == '__main__':
    count = {}
    reader = 'LLAMA'
    filenames = [file for file in glob.glob(f'input/{reader}/*.xlsx')]
    if reader == 'GPT' or reader == 'LLAMA':
        FILES = FILES[1:]
    for f in FILES:
        out_name = f

        df = pd.read_excel(f'input/interactions/{reader}/{out_name}_reading_BioRECIPE.xlsx', index_col=None).fillna('nan').astype(str)
        df = df.applymap(lambda x: x.lower().strip() if isinstance(x, str) else x)
        #df['Regulator Name'] = df['Regulator Name'].str.replace('-', '').str.replace(' ', '').str.replace('_', '')

        merged_df = merge_duplicates(df,
                                     ["Regulator Name",
                                      "Regulator Type",
                                      "Regulator ID",
                                      "Regulated Name",
                                      "Regulated Type",
                                      "Regulated ID",
                                      "Sign"
        ])

        interaction_attr = len(merged_df)
        ele_id_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Regulator Database'] if x.lower() != 'be' and \
                                                                merged_df.loc[row, 'Regulator ID'] != 'nan']) + \
                      len([x for row in range(len(merged_df))
                           for x in merged_df.loc[row, 'Regulated Database'] if x.lower() != 'be' and \
                           merged_df.loc[row, 'Regulator ID'] != 'nan'])
        subtype_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Regulator Subtype'] if x.lower() != 'nan']) + \
                       len([x for row in range(len(merged_df))
                            for x in merged_df.loc[row, 'Regulated Subtype'] if x.lower() != 'nan'])
        sign_attr = len([x for x in merged_df['Sign'].to_list() if x.lower() != 'nan'])
        site_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Site'] if x.lower() != 'nan'])
        cell_line_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Cell Line'] if x.lower() != 'nan'])
        cell_type_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Cell Type'] if x.lower() != 'nan'])
        tissue_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Tissue Type'] if x.lower() != 'nan'])
        organism_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Organism'] if x.lower() != 'nan'])
        statements_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Statements'] if x.lower() != 'nan'])

        connection_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Connection Type'] if x.lower() not in ['none', 'nan']])
        connection_type_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Connection Type'] if x.lower() in ['d', 'direct', 'positive', 'true']])
        location_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Regulated Compartment'] if x.lower() not in ['none', 'nan']])
        mechanism_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Mechanism'] if x.lower() not in ['none', 'nan']])
        phos_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Mechanism'] if x.lower() == 'phosphorylation'])
        paper_attr = len([x for row in range(len(merged_df))
                                for x in merged_df.loc[row, 'Paper IDs'] if x.lower() not in ['nan', '', 'none']])
        if 'interaction' not in count: count['interaction'] = {}
        if 'ele_id' not in count: count['ele_id'] = {}
        if 'subtype' not in count: count['subtype'] = {}
        if 'sign' not in count: count['sign'] = {}

        if 'connection' not in count: count['connection'] = {}
        if 'connection_type' not in count: count['connection_type'] = {}
        if 'location' not in count: count['location'] = {}
        if 'mechanism' not in count: count['mechanism'] = {}
        if 'phosphorylation' not in count: count['phosphorylation'] = {}

        if 'site' not in count: count['site'] = {}
        if 'cell_line' not in count: count['cell_line'] = {}
        if 'cell_type' not in count: count['cell_type'] = {}
        if 'tissue' not in count: count['tissue'] = {}
        if 'organism' not in count: count['organism'] = {}
        if 'statements' not in count: count['statements'] = {}
        if 'paper_id' not in count: count['paper_id'] = {}

        count['interaction'][out_name] = interaction_attr
        count['ele_id'][out_name] = ele_id_attr
        count['subtype'][out_name] = subtype_attr
        count['sign'][out_name] = sign_attr
        count['site'][out_name] = site_attr
        count['cell_line'][out_name] = cell_line_attr
        count['cell_type'][out_name] = cell_type_attr
        count['tissue'][out_name] = tissue_attr
        count['organism'][out_name] = organism_attr
        count['statements'][out_name] = statements_attr
        count['connection'][out_name] = connection_attr
        count['connection_type'][out_name] = connection_type_attr
        count['location'][out_name] = location_attr
        count['mechanism'][out_name] = mechanism_attr
        count['phosphorylation'][out_name] = phos_attr
        count['paper_id'][out_name] = paper_attr

    dict_ = {}
    dict_['reading'] = FILES
    for key, value in count.items():
        dict_[key] = list(count[key].values())

    count_df = pd.DataFrame(dict_)

    count_df.to_csv(f'{reader}_summary.csv', index=False)


