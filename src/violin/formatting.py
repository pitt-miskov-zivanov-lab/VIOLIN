"""
formatting.py

Handles the model and reading formatting functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
Updated May 2024 - Haomiao Luo
"""

import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import os.path
import logging
import re

# define regex for valid characters in variable names
_VALID_CHARS = r'a-zA-Z0-9\_'

# valid element types
_VALID_TYPES = [
    'protein', 'protein family', 'protein complex',
    'rna', 'mrna', 'gene', 'chemical', 'biological process'
    ]

_VAR_COL = 'Variable'
_IDX_COL = '#'
required_model = ['Element Name','Element Type','Element IDs','Variable','Positive Regulator List', 'Negative Regulator List']

type_abbr_dict = {
             'proteinfamily': 'pf',
             'proteincomplex': 'pf',
             'protein': 'pn',
             'chemical': 'che',
             'chemicalfamily': 'cf',
             'gene': 'gene',
             'rna': 'rna',
             'mutation': 'mut',
             'biologicalprocess': 'bp'
}

subtype_abbr_dict = {
                "enzyme":"enz",
                "transporter":"tsp",
                "transcription-factor":"tsf",
                "transcription-repressor":"tsr",
                "transducer":"tsd",
                "kinase":"kin",
                "interferon":"ifn",
                "interleukin":"ilk",
                "subunit":"sub",
                "cytokine":"cyt",
                "tyrosine":"tyr",
                "receptor":"rec",
                "caspase":"cas",
                "phosphatase":"pho",
                "adaptor":'ada',
                "peptidase":"pep",
                "cyclin":"cyc",
                "growth-factor":"gwf",
                "binding":"bin",
                "molecule":"mol",
                "oncogene":"onc",
                "proto-oncogene":"pnc",
                "suppressor":"sup",
                "tumor":"tum",
                "signaling":"sig",
                "biological":"bio",
                "process":"prc",
                "protein":"prt",
                "redox":"red",
                "metallopeptidase":"mtp",
                "nonhistone":"nhs",
                "nucleoprotein":"ncp",
                "hormone":"hor",
                "ligase":'lgs',
                "ligand":"lgd",
                "regulator":"reg",
                "ubiquitin-protein":"ubp",
                "catalytic":"cat",
                "gtpase":"gtp",
                "reverse":"rvs",
                "transcriptase":"tst",
                "dehydrogenase":"dhy",
                "hydrogenase":"hyd",
                "peroxidase":"pox",
                "oxidase":"oxi",
                "glycoprotein":"glp",
                "necrosis-factor":"nec",
                "apoptosis":"apo",
                "active":"act"
}


def evidence_score(reading_df, col_names):
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

    #Counting the number of duplicates
    counted_reading['Evidence Score'] = counted_reading[remainder[0]].str.len()
    # counted_reading.to_csv("Trying.csv")

    return counted_reading

def add_regulator_names_id(model_df):
    """
    This function converts the model regulator lists from BioRECIPE variables to the common element names and database identifiers

    Parameters
    ----------
    model_df : pd.DataFrame
        The model dataframe (in BioRECIPE format)

    Returns
    -------
    model_df : pd.DataFrame
        A new dataframe with added columns containing the positive and negative regulators listed
        by their Element Names and IDs
    """
    #removes the initial values from the model dataframe, as they're not needed
    #Also adds new columns for the positive and negative regulator names and IDs
    col_headers = list(model_df.columns)
    model_df = model_df[col_headers]
    reg_col_list = ['Positive Regulator List', 'Negative Regulator List']
    model_df[reg_col_list] = model_df[reg_col_list].apply(lambda x: x.astype(str).str.lower())
    #Columns for positive
    model_df['Positive Names'] = pd.Series().astype(object)
    model_df['Positive IDs'] = pd.Series().astype(object)
    model_df['Negative Names'] = pd.Series().astype(object)
    model_df['Negative IDs'] = pd.Series().astype(object)

    #Convert Regulators
    for sign in ['Negative','Positive']:
        for y in range(model_df.shape[0]):
            if model_df[sign+' Regulator List'][y] in ['', "nan"]:
                model_df.at[y,sign+' Names'] = "nan"
                model_df.at[y,sign+' IDs'] = "nan"

            else:
                reg_name = model_df[sign+" Regulator List"][y].split(",")
                if '' in reg_name: reg_name.remove('')
                reg_id = []
                reg_var = reg_name.copy()
                model_df.at[y,sign+' Regulator List'] = reg_var

                #find index for regulator in variable column, and copy the Element Name and IDs to the new columns
                for element in reg_name:
                    idx = list(model_df['Listname']).index(element)
                    reg_name[reg_name.index(element)] = model_df['Element Name'][idx]
                    #idx = list(model_df["Element Name"]).index(element)
                    #Since there are multiple IDs for each element, need to keep track of which
                    #IDs go with which regulator
                    reg_id.append(model_df["Element IDs"][idx])

                model_df.at[y,sign+' Names'] = reg_name
                model_df.at[y,sign+' IDs'] = reg_id

    return model_df
def format_variable_names(model: pd.DataFrame):

    """
    This function formats model variable names to make compatible with model checking

    Parameters
    ----------
    model: DataFrame
        A dataframe of model file
    Returns
    -------
    model: DataFrame
        model dataframe with standardized variable names

    """

    global _VALID_CHARS
    global _VAR_COL

    # remove whitespace in variable names
    model[_VAR_COL] = model[_VAR_COL].str.strip()

    # collect invalid element names in a list so they can be removed everywhere in the model
    # find invalid characters in element names and names starting with numbers
    invalid_names = [
        x for x in model[_VAR_COL]
        if re.search(r'(^[0-9]+)',x.strip()) or re.search(r'([^'+_VALID_CHARS+']+)',x.strip())
        ]

    if len(invalid_names) > 0:
        logging.info('Formatting variable names: ')

    # remove invalid characters at the start of the variable name
    replace_names = [re.sub(r'^[^'+_VALID_CHARS+']+','',x) for x in invalid_names]
    # replace invalid characters elsewhere in variable names
    replace_names = [re.sub(r'[^'+_VALID_CHARS+']+','_',x) for x in replace_names]

    # add ELE_ at the beginning of names starting with numbers
    replace_names = [re.sub(r'(^[0-9]+)','ELE_\\1',x) for x in replace_names]

    name_pairs = zip(invalid_names,replace_names)

    for (invalid_name,replace_name) in name_pairs:
        logging.info('%s -> %s' % (invalid_name,replace_name))
        model.replace(re.escape(invalid_name),re.escape(replace_name),regex=True,inplace=True)

    return model

def get_type(input_type):

    """
    This function standardizes element types

    Parameters
    ----------
    input_type: str
        A string of entity type
    Returns
    -------
    standard string to describe the type of entity

    """

    global _VALID_TYPES

    if input_type.lower() in _VALID_TYPES:
        return input_type
    elif input_type.lower().startswith('protein'):
        return 'protein'
    elif input_type.lower().startswith('chemical'):
        return 'chemical'
    elif input_type.lower().startswith('biological'):
        return 'biological'
    else:
        return 'other'


def get_element(reg_rule, layer):

    """
    This function parses the regulation rule and disentangle the symbol operators converting rule to a list of regulators

    Parameters
    ----------
    reg_rule: str
        A BioRECIPE Regulation Rule
    layer: str
        counter for recursive time, the default is 0
    Returns
    -------
    regulator_list: list
        A list of regulators

    """

    if reg_rule:
        regulator_list = []

        if '+' not in reg_rule:
            reg_list = split_comma_out_parentheses(reg_rule)
        else:
            if ',' in reg_rule:
                raise ValueError(
                'Found mixed commas and plus sign in regulation function'
                )
            elif reg_rule[-1] == '+':
                raise ValueError(
                'Regulation rule is not correct'
                )
            else:
                reg_list = reg_rule.split('+')

        for reg_element in reg_list:
            if reg_element[0] == '{' and reg_element[-1] == '}':
                assert(layer == 0)
                if '*' in reg_element:
                    weight, name = reg_element[1:-1].split('*')
                    regulator_list = regulator_list + get_element(name, 1)
                else:
                    regulator_list = regulator_list + get_element(reg_element, 1)

            elif reg_element[0] == '{' and reg_element[-1] == ']':
                # This is a necessary pair
                # check the point between {} and []
                parentheses = 0
                cutpoint = 0
                for index, char in enumerate(reg_element):
                    if char == '{':
                        parentheses +=1
                    elif char =='}':
                        parentheses -=1

                    if parentheses == 0:
                        cutpoint = index
                        break

                necessary_element = reg_element[1: cutpoint]
                enhence_element = reg_element[cutpoint+2:-1]

                if '*' in necessary_element:
                    weight, name = necessary_element.split('*')
                    regulator_list = regulator_list + get_element(name, 1)
                else:
                    regulator_list = regulator_list + get_element(necessary_element, 1)

                if '*' in enhence_element:
                    weight, name = enhence_element.split('*')
                    regulator_list = regulator_list + get_element(name, 1)
                else:
                    regulator_list = regulator_list + get_element(enhence_element, 1)

            elif reg_element[0] == '(' and reg_element[-1] == ')':
                list = [element for ele_list in split_comma_out_parentheses(reg_element[1:-1])
                        for element in get_element(ele_list, 1)]
                regulator_list += list
            else:
                #print(f'reg_element: {reg_element}')
                assert(',' not in reg_element)

                if reg_element[-1] == '^':
                    regulator_list.append(reg_element[0:-1])
                elif '&' in reg_element:
                    regulator_list.append(reg_element[1:-1])
                elif '*' in reg_element:
                    multiply_reg_list = reg_element.split('*')
                    for reg_ in multiply_reg_list:
                        if re.search(r'^[0-9]', reg_):
                            pass
                        elif not re.search(r'[a-zA-Z0-9\_!]+', reg_):
                            pass
                        else:
                            regulator_list.append(reg_)
                elif reg_element[0] == '!':
                    if '~' in reg_element[1:]:
                        delay, reg_delay = reg_element[1:].split('~')
                        regulator_list.append(reg_delay)
                    else:
                        regulator_list.append(reg_element[1:])

                elif '=' in reg_element:
                    name, target_state = reg_element.split('=')
                    regulator_list.append(target_state)
                elif '~' in reg_element:
                    delay, state = reg_element.split('~')
                    regulator_list.append(state)

                else:
                    regulator_list.append(reg_element)

        return regulator_list


def split_comma_out_parentheses(reg_rule):

    """
    This function split the parentheses by comma outside of parentheses. e.g. '(A,B),(C,B)' -> ['(A,B)','(C,B)']

    Parameters
    ----------
    reg_rule: str
        A regulation rule

    Returns
    -------
    reg_list: list
    A list of expressions that are separated by brackets

    """

    reg_list = list()
    parentheses = 0
    start = 0
    for index, char in enumerate(reg_rule):
        if index == len(reg_rule) - 1:
            reg_list.append(reg_rule[start:index+1])
        elif char == '(' or char == '{' or char == '[':
            parentheses += 1
        elif char == ')' or char == '}' or char == ']':
            parentheses -= 1
        elif (char == ',' and parentheses == 0):
            reg_list.append(reg_rule[start:index])
            start = index+1
    return reg_list


def wrap_list_to_str(df, cols):
    """
    This function wraps the lists in the output dataframe to strings
    Parameters
    ----------
    df: pd.DataFrame
        output dataframe
    cols: list
        a list of columns name
    Returns
    -------
    df: pd.DataFrame
    """
    for row in range(len(df)):
        for col in cols:
            df.loc[row, col] = ','.join(list(df.loc[row, col]))
    return df


def get_listname(idx, model_df):
    """
        Create the list-names by element attributes
    Parameters
    ----------
    idx: int
        the index of element
    model_df: pd.DataFrame
        the model table
    Returns
    -------
    listname: str
        formatted name for regulator list column
    """
    ele_col_list = ['Element Name', 'Element Type', 'Element Subtype', 'Compartment ID']
    model_df[ele_col_list] = model_df[ele_col_list].apply(lambda x: x.astype(str).str.lower())
    if str(model_df.loc[idx, 'Element Type']).replace(' ', '') not in type_abbr_dict:
        ele_type = model_df.loc[idx, 'Element Type'].replace(' ', '')
    else:
        ele_type = type_abbr_dict[model_df.loc[idx, 'Element Type'].replace(' ', '')]
    listname = '{}_{}_{}_{}'.format(
        model_df.loc[idx, 'Element Name'],
        ele_type,
        get_subtype_abbr(model_df.loc[idx, 'Element Subtype']),
        model_df.loc[idx, 'Compartment ID'].replace(':', '')
    )
    return listname


def get_subtype_abbr(subtype):
    """

    Parameters
    ----------
    subtype: str
        The subtype of the element
    Returns
    -------
    abbr: str
        The abbreviation of the first subtype
    """
    list_ = []
    # Only get first subtype (TBD for the other subtypes)
    subtype = subtype.split(',')[0]
    if subtype not in ['', 'nan']:
        for x in [subname for subname in subtype.replace('(', ' ').replace(')', ' ').split(' ') if
                  subname not in ['', ' ']]:
            if x.lower() not in subtype_abbr_dict.keys():
                list_.append(x.strip())
            else:
                list_.append(subtype_abbr_dict[x.lower().strip()])
        abbr = ''.join(list_)
    else:
        abbr = 'nan'
    return abbr

#TODO: implement with functionality and integrate with BioRECIPE
# def convert_to_biorecipes(model, att_list = [], separate = True):
#     """
#     This function imports a model which is NOT in the BioRECIPES format,
#     such as models formatted as node-edge lists.
#     Regulators may be represented in the REACH formatt, separated by regulator sign,
#     or unseparated, with a speicifed column for regulator sign
#
#     Parameters
#     ----------
#     model : str
#         Directory and filename of the file containing the model BioRECIPES spreadsheet
#         Accepted files: .txt, .csv, .tsv, .xlsx
#     model_cols : list
#         Column names of the model file. Default names are found in required_model
#     att_list : list
#         List of Element attributes (in addition to Name, ID, and Type)
#         Default is no additional attributes
#     separate : Boolean
#         Whether or not the model presents regulator in separate Positive/Negative columns (True)
#         or in a single column with Regulator Sign attribute (False)
#         Default is True
#
#     Returns
#     -------
#     new_model : pd.DataFrame
#         Formatted model dataframe
#     """
#     #Upload the model files as dataframes based on the file extension
#     model_ext = os.path.splitext(model)[1]
#
#     if model_ext == '.txt': model_df = pd.read_csv(model, sep='\t',index_col=None).fillna("nan")
#     elif model_ext == '.csv': model_df = pd.read_csv(model, sep=',',index_col=None).fillna("nan")
#     elif model_ext == '.xlsx': model_df = pd.read_excel(model, index_col=None).fillna("nan")
#     elif model_ext == '.tsv': model_df = pd.read_csv(model, sep='\t',index_col=None).fillna("nan")
#     else: raise ValueError("The accepted file extensions are .txt, .csv, .xlsx, and .tsv")
#
#     #Get column header names as list
#     col_names = list(model_df.columns)
#
#     #formate to bare minimum information
#     bare_cols = [x.lower().replace(" ","").replace("_","").replace("-","") for x in col_names]
#
#     #Accepted target/regulated headers
#     t_name_list = ["elementname","targetname","regulatedname"]
#     t_type_list = ["elementtype","targettype","regulatedtype"]
#     t_id_list = ["elementid","elementidentifier","targetid","targetidentifier","regulatedid","regulatredidentifier"]
#     t_att_pre = ["element","target","regulated",""]
#
#     #If the variables are separated into Positive/Negative:
#     if separate:
#         #Accepted source/regulator headers (assuming separate positive and negative columns)
#         s_name_list = ["positiveregname","posregname","positiveregulatorname","posregulatorname",
#                         "positivesourcename","possourcename"]
#         s_type_list = ["positiveregtype","posregtype","positiveregulatortype","positiveregtype",
#                         "positivesourcetype","possourcetype"]
#         s_id_list = ["positiveregid","posregid","positiveregulatorid","positiveregid",
#                         "positivesourceid","possourceid","positiveregidentifier","posregidentifier",
#                         "positiveregulatoridentifier","positiveregidentifier","positivesourceidentifier","possourceidentifier"]
#         s_att_pre = ["posreg","posregulator","positivereg","positiveregulator"]
#
#         #Check intersection of accepted column names and file column names
#         if {len(set(t_name_list) & set(bare_cols)) == 1 & len(set(t_type_list) & set(bare_cols)) == 1 &
#             len(set(t_id_list) & set(bare_cols)) == 1 & len(set(s_name_list) & set(bare_cols)) == 1 &
#             len(set(s_type_list) & set(bare_cols)) == 1 & len(set(s_id_list) & set(bare_cols)) == 1}:
#             #If minimum necessary columns are found, define variables for the column header
#             target_name = col_names[bare_cols.index((set(t_name_list) & set(bare_cols)).pop())]
#             target_type = col_names[bare_cols.index((set(t_type_list) & set(bare_cols)).pop())]
#             target_id = col_names[bare_cols.index((set(t_id_list) & set(bare_cols)).pop())]
#             pos_source_name = col_names[bare_cols.index((set(s_name_list) & set(bare_cols)).pop())]
#             pos_source_type = col_names[bare_cols.index((set(s_type_list) & set(bare_cols)).pop())]
#             pos_source_id = col_names[bare_cols.index((set(s_id_list) & set(bare_cols)).pop())]
#             neg_source_name = pos_source_name.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
#             neg_source_type = pos_source_type.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
#             neg_source_id = pos_source_id.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
#             #store column header names in a dictionary
#             model_cols = {"target_name" : target_name,"target_type" :target_type,"target_id" : target_id,
#                             "pos_source_name" : pos_source_name, "pos_source_type" : pos_source_type, "pos_source_id" : pos_source_id,
#                             "neg_source_name" : neg_source_name, "neg_source_type" : neg_source_type, "neg_source_id" : neg_source_id}
#             # Add Connection Type when available
#             if "Connection Type" in col_names: model_cols['cxn_type'] = "Connection Type"
#             #Now for the attributes:
#             for x in att_list:
#                 t_att_list = [pre + x.lower().replace(" ","") for pre in t_att_pre]
#                 s_att_list = [pre + x.lower().replace(" ","") for pre in s_att_pre]
#                 #made sure the attribute is in the reading columns
#                 if len(set(t_att_list) & set(bare_cols)) == 1 & len(set(s_att_list) & set(bare_cols)) == 1:
#                     #if it's found, add it to the reading columns
#                     #add the attribute for the target
#                     model_cols['target_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(t_att_list) & set(bare_cols)).pop())]
#                     model_cols['pos_source_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(s_att_list) & set(bare_cols)).pop())]
#                     model_cols['neg_source_'+x.lower().replace(" ","_")] = model_cols['pos_source_'+x.lower().replace(" ","_")].replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
#                 else:
#                     raise ValueError("Attribute \""+x+"\" was not found in your Model input document."+"\n"+
#                     "Please check your file and try again")
#             #Import the Element and regulator sets
#             #Both regulator sets will have "nan" items, representing those rows which do not have
#             #a regulator of that sign
#             elements = set(model_df[model_cols['target_name']])
#             pos_regs = set(model_df[model_cols['pos_source_name']])
#             pos_regs.remove('nan')
#             neg_regs = set(model_df[model_cols['neg_source_name']])
#             neg_regs.remove('nan')
#
#             #find regulators which are not already in the element list
#             pos_not_elements = np.setdiff1d(pos_regs,elements)[0]
#             neg_not_elements = np.setdiff1d(neg_regs,elements)[0]
#
#             #get column names of all regulator attributes
#             pos_cols = {key: value for key, value in model_cols.items() if 'pos' in key.lower()}
#             neg_cols = {key: value for key, value in model_cols.items() if 'neg' in key.lower()}
#
#             #get indices of not_elements
#             pos_idx = []
#             neg_idx = []
#             for x in list(pos_not_elements):
#                 pos_idx.append(list(model_df[model_cols['pos_source_name']]).index(x))
#             for y in list(neg_not_elements):
#                 neg_idx.append(list(model_df[model_cols['neg_source_name']]).index(y))
#
#             #create subsets of DF based on not_elements
#             pos_sub = model_df.loc[pos_idx,pos_cols.values()]
#             neg_sub = model_df.loc[neg_idx,neg_cols.values()]
#
#             #Reduce sub DF to only unique rows
#             unique_pos = pos_sub.drop_duplicates()
#             unique_neg = neg_sub.drop_duplicates()
#
#             #Add subsets to model_df
#             for each in list(pos_cols.keys()):
#                 unique_pos.rename(columns={pos_cols[each]: model_cols[each.replace('pos_source','target')]},inplace=True)
#             for every in list(neg_cols.keys()):
#                 unique_neg.rename(columns={neg_cols[every]: model_cols[every.replace('neg_source','target')]},inplace=True)
#             model_df = model_df.append(unique_pos,ignore_index=True).append(unique_neg,ignore_index=True).fillna('nan')
#
#
#             # Change column header names
#             biorecipes_cols = {'target_name':'Element Name', 'target_type':'Element Type',
#                                'target_id':'Element IDs','pos_source_name':'Positive Regulator List',
#                                'neg_source_name':'Negative Regulator List'}
#             for x in list(biorecipes_cols.keys()):
#                 model_df = model_df.rename(columns={model_cols[x]:biorecipes_cols[x]})
#
#             #Delete Extraneous Columns
#             model_df = model_df.drop(columns=list(set(pos_cols.values())&set(model_df.columns)))
#             model_df = model_df.drop(columns=list(set(neg_cols.values())&set(model_df.columns)))
#
#             group_cols = [value for key, value in biorecipes_cols.items() if 'target' in key.lower()]
#             remainder = [x for x in model_df.columns if x not in group_cols]
#
#             #As VIOLIN Identifies duplicates, it merges attributes from the remainder list into a single cell
#             biorecipes_model = model_df.groupby(group_cols)[remainder[0]].apply(list).reset_index(name=remainder[0])
#             for x in range(1,len(remainder)):
#                 sub = model_df.groupby(group_cols)[remainder[x]].apply(list).reset_index(name=remainder[x])
#                 biorecipes_model[remainder[x]] = sub[remainder[x]]
#             for each in remainder:
#                 biorecipes_model[each] = biorecipes_model[each].apply(','.join)
#             biorecipes_model = biorecipes_model.replace({',nan': ''}, regex=True)
#             biorecipes_model = biorecipes_model.replace({'nan': ''}, regex=True)
#
#             biorecipes_model = biorecipes_model.sort_values(by='Element Name',ascending=True)
#
#
#             #If Variables present, covert regulator variable name lists to common names and database identifiers
#             if 'Variable' not in list(model_cols.keys()):
#                 biorecipes_model['Variable'] = biorecipes_model['Element Name']
#                 # remove whitespace from variable names
#                 biorecipes_model['Variable'] = biorecipes_model['Variable'].map(lambda x: x.lstrip(' '))
#                 biorecipes_model = add_regulator_names_id(biorecipes_model)
#             else: biorecipes_model = add_regulator_names_id(biorecipes_model)
#         else:
#             raise ValueError("Unaccepted Column Names. Please check that you have"+"\n"+
#                              "Names, Types, and IDs for both source and target nodes")
#
#     #If variables are presented in a single column, with a "regulator sign" column
#     else:
#                 #Accepted source/regulator headers (assuming separate positive and negative columns)
#         s_name_list = ["regname","regulatorname","sourcename"]
#         s_type_list = ["regtype","regulatortype","sourcetype"]
#         s_id_list = ["regid","regulatorid","sourceid",
#                      "regidentifier","regulatoridentifier","sourceidentifier"]
#         s_sign_list = ["regsign","regulatorsign","regulationsign"]
#         s_att_pre = ["reg","regulator","source"]
#
#         #Check intersection of accepted column names and file column names
#         if {len(set(t_name_list) & set(bare_cols)) == 1 & len(set(t_type_list) & set(bare_cols)) == 1 &
#             len(set(t_id_list) & set(bare_cols)) == 1 & len(set(s_name_list) & set(bare_cols)) == 1 &
#             len(set(s_type_list) & set(bare_cols)) == 1 & len(set(s_id_list) & set(bare_cols)) == 1 &
#             len(set(s_sign_list) & set(bare_cols)) == 1}:
#             #If minimum necessary columns are found, define variables for the column header
#             target_name = col_names[bare_cols.index((set(t_name_list) & set(bare_cols)).pop())]
#             target_type = col_names[bare_cols.index((set(t_type_list) & set(bare_cols)).pop())]
#             target_id = col_names[bare_cols.index((set(t_id_list) & set(bare_cols)).pop())]
#             source_name = col_names[bare_cols.index((set(s_name_list) & set(bare_cols)).pop())]
#             source_type = col_names[bare_cols.index((set(s_type_list) & set(bare_cols)).pop())]
#             source_id = col_names[bare_cols.index((set(s_id_list) & set(bare_cols)).pop())]
#             source_sign = col_names[bare_cols.index((set(s_sign_list) & set(bare_cols)).pop())]
#             #store column header names in a dictionary
#             model_cols = {"target_name" : target_name,"target_type" :target_type,"target_id" : target_id,
#                             "source_name" : source_name, "source_type" : source_type, "source_id" : source_id,
#                             "regulation_sign" : source_sign}
#             # Add Connection Type when available
#             if "Connection Type" in col_names: model_cols['cxn_type'] = "Connection Type"
#
#             #Now for the attributes:
#             for x in att_list:
#                 t_att_list = [pre + x.lower().replace(" ","") for pre in t_att_pre]
#                 s_att_list = [pre + x.lower().replace(" ","") for pre in s_att_pre]
#                 #made sure the attribute is in the reading columns
#                 if len(set(t_att_list) & set(bare_cols)) == 1 & len(set(s_att_list) & set(bare_cols)) == 1:
#                     #if it's found, add it to the reading columns
#                     #add the attribute for the target
#                     model_cols['target_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(t_att_list) & set(bare_cols)).pop())]
#                     model_cols['source_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(s_att_list) & set(bare_cols)).pop())]
#                 else:
#                     raise ValueError("Attribute \""+x+"\" was not found in your LEE input document."+"\n"+
#                     "Please check your file and try again")
#
#             #Import the Element and regulator sets
#             elements = set(model_df[model_cols['target_name']])
#             regs = set(model_df[model_cols['source_name']])
#             regs.remove('nan')
#
#             #find regulators which are not elements
#             not_elements = np.setdiff1d(regs,elements)[0]
#
#             #get column names of all regulator attributes
#             reg_col_names = [value for key, value in model_cols.items() if 'source' in key.lower()]
#
#             #get indices of not_elements
#             not_idx = []
#             for x in list(not_elements):
#                 not_idx.append(list(model_df[model_cols['source_name']]).index(x))
#
#             #create subsets of DF based on not_elements
#             reg_sub = model_df.loc[not_idx,reg_col_names]
#
#             #Reduce sub DF to only unique rows
#             unique_reg = reg_sub.drop_duplicates()
#
#             #Add subsets to model_df
#             unique_reg.columns = [s.replace(model_cols['source_name'], model_cols['target_name']) for s in unique_reg.columns]
#             model_df = model_df.append(unique_reg,ignore_index=True).fillna('nan')
#
#             #Need to split the regulators into positive/negative before combining rows
#             if "Connection Type" in col_names: reg_atts = reg_col_names + [model_cols['cxn_type']]
#             else: reg_atts = reg_col_names
#
#             for x in range(len(list(model_df[model_cols['source_name']]))):
#                 if model_df.loc[x,model_cols['regulation_sign']].lower() in ['negative','decrease','inhibit']:
#                     for each in reg_atts:
#                         model_df.loc[x,'Negative '+each] =  model_df.loc[x,each]
#                 else:
#                     for each in reg_atts:
#                         model_df.loc[x,'Positive '+each] =  model_df.loc[x,each]
#             model_df = model_df.drop(columns=reg_atts+[model_cols['regulation_sign']]).fillna('nan')
#
#             #Now need to combine rows with the same element
#             #remainders are the regulator columns which need to be retained
#             if "Connection Type" in col_names: remainder = ['Positive ' + s for s in reg_col_names] + ['Positive Connection Type List']+['Negative ' + s for s in reg_col_names]+['Negative Connection Type List']
#             else: remainder = ['Positive ' + s for s in reg_col_names] +['Negative ' + s for s in reg_col_names]
#             model_cols = [x for x in model_df if x not in remainder]
#
#             # Change column header names
#             biorecipes_cols = {'target_name':'Element Name', 'target_type':'Element Type',
#                                'target_id':'Element IDs'}
#             for x in biorecipes_cols:
#                 model_df = model_df.rename(columns={model_cols[x]:biorecipes_cols[x]})
#
#             #As VIOLIN Identifies duplicates, it merges attributes from the remainder list into a single cell
#             biorecipes_model = model_df.groupby(model_cols)[remainder[0]].apply(list).reset_index(name=remainder[0])
#             for x in range(1,len(remainder)):
#                 sub = model_df.groupby(model_cols)[remainder[x]].apply(list).reset_index(name=remainder[x])
#                 biorecipes_model[remainder[x]] = sub[remainder[x]]
#             for each in remainder:
#                 biorecipes_model[each] = biorecipes_model[each].apply(','.join)
#             biorecipes_model = biorecipes_model.replace({',nan': ''}, regex=True)
#             biorecipes_model = biorecipes_model.replace({'nan': ''}, regex=True)
#
#             biorecipes_model = biorecipes_model.sort_values(by='Element Name',ascending=True)
#
#             #If Variables present, covert regulator variable name lists to common names and database identifiers
#             if 'Variable' not in model_cols:
#                 biorecipes_model['Variable'] = biorecipes_model['Element Name']+'_'+biorecipes_model['Element Type']
#                 # remove whitespace from variable names
#                 biorecipes_model['Variable'] = biorecipes_model['Variable'].map(lambda x: x.lstrip(' '))
#                 biorecipes_model = add_regulator_names_id(biorecipes_model)
#             else: biorecipes_model = add_regulator_names_id(biorecipes_model)
#         else:
#             raise ValueError("Unaccepted Column Names. Please check that you have"+"\n"+
#                              "Names, Types, and IDs for both source and target nodes")
#
#         biorecipes_model = biorecipes_model.apply(lambda x: x.astype(str).str.lower())
#
#     return biorecipes_model

#TODO: implement the functionality and integrate with BioRECIPE
# def convert_reading(reading, action, atts = []):
#     """
#     This function formats the machine reading output,
#     either separating regulator names and attributes into 'positive' and 'negative' columns to match REACH formatting,
#     or combining regulator names and attributes without regulator sign distinction, and adding a 'regulator sign' column.
#     This function can take the machine reading as either a filename or as an already uploaded dataframe.
#
#     Parameters
#     ----------
#     reading : str or pd.DataFrame
#         Machine reading output, either as file location string or dataframe
#     action : str
#         Action to be performed by function
#         Accepts only 'combine' or 'separate' as input
#     atts : list
#         List of attributes associated with each regualtor
#         Default list is ['Type','ID']
#         List should not include regulator signs (where applicable)
#
#     Returns
#     -------
#     reading_df : pd.DataFrame
#         A dataframe with the specified formatting completed
#     """
#
#
#     if type(reading) == str:
#         reading_ext = os.path.splitext(reading)[1]
#         if reading_ext == '.txt': reading_df = pd.read_csv(reading, sep='\t',index_col=None).fillna("nan")
#         elif reading_ext == '.csv': reading_df = pd.read_csv(reading, sep=',',index_col=None).fillna("nan")
#         elif reading_ext == '.xlsx': reading_df = pd.read_excel(reading,index_col=None).fillna("nan")
#         elif reading_ext == '.tsv': reading_df = pd.read_csv(reading, sep='\t',index_col=None).fillna("nan")
#         else: raise ValueError("The accepted file extensions are .txt, .csv, .xslx, and .tsv")
#
#     elif type(reading) == pd.DataFrame:
#         reading_df = reading
#
#     else: raise ValueError("Unsupported input type. This functions accepts filenames and dataframes")
#
#
#     if action == 'separate':
#         #Begin relative column name retrieval
#         #Accepted source/regulator headers (assuming separate positive and negative columns)
#         s_name_list = ["regname","regulatorname","sourcename"]
#         s_type_list = ["regtype","regulatortype","sourcetype"]
#         s_id_list = ["regid","regulatorid","sourceid","regidentifier",
#                         "regulatoridentifier","sourceidentifier"]
#         s_att_pre = ["reg","regulator","source"]
#         s_sign_list = "sign","regsign","regulatorsign","regulationsign"
#
#         #Get column header names as list
#         col_names = list(reading_df.columns)
#
#         #formate to bare minimum information
#         bare_cols = [x.lower().replace(" ","").replace("_","").replace("-","") for x in col_names]
#
#         #Check intersection of accepted column names and file column names
#         if {len(set(s_name_list) & set(bare_cols)) == 1 &
#             len(set(s_type_list) & set(bare_cols)) == 1 & len(set(s_id_list) & set(bare_cols)) == 1}:
#             #If minimum necessary columns are found, define variables for the column header
#             source_name = col_names[bare_cols.index((set(s_name_list) & set(bare_cols)).pop())]
#             source_type = col_names[bare_cols.index((set(s_type_list) & set(bare_cols)).pop())]
#             source_id = col_names[bare_cols.index((set(s_id_list) & set(bare_cols)).pop())]
#             source_sign = col_names[bare_cols.index((set(s_sign_list) & set(bare_cols)).pop())]
#             #store column header names in a dictionary
#             reading_cols = {"source_name" : source_name, "source_type" : source_type,
#                             "source_id" : source_id, "source_sign" : source_sign}
#             #Now for the attributes:
#             for x in atts:
#                 s_att_list = [pre + x.lower().replace(" ","") for pre in s_att_pre]
#                 #made sure the attribute is in the reading columns
#                 if len(set(s_att_list) & set(bare_cols)) == 1:
#                     #if it's found, add it to the reading columns
#                     #add the attribute for the target
#                     reading_cols['source_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(s_att_list) & set(bare_cols)).pop())]
#                 else:
#                     raise ValueError("Attribute \""+x+"\" was not found in your LEE input document."+"\n"+
#                     "Please check your file and try again")
#         else:
#             raise ValueError("Your LEE input is missing information."+"\n"+
#             "VIOLIN requires the following information: Name, Type, and ID of target node and regulators")
#         #End relative column name retrieval
#
#         #make sure necessary header column headers are present
#
#         #Need to split the regulators into positive/negative before combining rows
#         new_reading_cols = {}
#         for x in range(len(list(reading_df[reading_cols["source_name"]]))):
#             if reading_df.loc[x,reading_cols["source_sign"]].lower() in ['negative','decrease','decreases','inhibit']:
#                 for each in list(reading_cols.keys()):
#                     reading_df.loc[x,'Negative '+reading_cols[each]] =  reading_df.loc[x,reading_cols[each]]
#                     new_reading_cols['neg_'+each] = 'Negative '+reading_cols[each]
#             else:
#                 for each in list(reading_cols.keys()):
#                     reading_df.loc[x,'Positive '+reading_cols[each]] =  reading_df.loc[x,reading_cols[each]]
#                     new_reading_cols['pos_'+each] = 'Positive '+reading_cols[each]
#         #Delete regulator sign column and unsigned columns
#         for key in reading_cols:
#                 reading_df = reading_df.drop(columns=reading_cols[key]).fillna('nan')
#
#         #Add target names to new_reading_cols
#         #Accepted target/regulated headers
#         t_name_list = ["elementname","targetname","regulatedname"]
#         t_type_list = ["elementtype","targettype","regulatedtype"]
#         t_id_list = ["elementid","elementidentifier","targetid","targetidentifier","regulatedid","regulatredidentifier"]
#         t_att_pre = ["element","target","regulated",""]
#         if {len(set(t_name_list) & set(bare_cols)) == 1 & len(set(t_type_list) & set(bare_cols)) == 1 &
#         len(set(t_id_list) & set(bare_cols)) == 1}:
#             #If minimum necessary columns are found, define variables for the column header
#             target_name = col_names[bare_cols.index((set(t_name_list) & set(bare_cols)).pop())]
#             target_type = col_names[bare_cols.index((set(t_type_list) & set(bare_cols)).pop())]
#             target_id = col_names[bare_cols.index((set(t_id_list) & set(bare_cols)).pop())]
#             new_reading_cols['target_name'] = target_name
#             new_reading_cols['target_type'] = target_type
#             new_reading_cols['target_id'] = target_id
#             #And the attributes
#             for x in atts:
#                 t_att_list = [pre + x.lower().replace(" ","") for pre in t_att_pre]
#                 if len(set(t_att_list) & set(bare_cols)) == 1:
#                     reading_cols['target_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(t_att_list) & set(bare_cols)).pop())]
#                 else:
#                     raise ValueError("Attribute \""+x+"\" was not found in your LEE input document."+"\n"+
#                     "Please check your file and try again")
#
#         return reading_df, new_reading_cols
#
#
#     elif action == 'combine':
#         #Begin relative column name retrieval
#         #Accepted source/regulator headers (assuming separate positive and negative columns)
#         s_name_list = ["positiveregname","posregname","positiveregulatorname","posregulatorname",
#                         "positivesourcename","possourcename"]
#         s_type_list = ["positiveregtype","posregtype","positiveregulatortype","positiveregtype",
#                         "positivesourcetype","possourcetype"]
#         s_id_list = ["positiveregid","posregid","positiveregulatorid","positiveregid",
#                         "positivesourceid","possourceid","positiveregidentifier","posregidentifier",
#                         "positiveregulatoridentifier","positiveregidentifier","positivesourceidentifier","possourceidentifier"]
#         s_att_pre = ["posreg","posregulator","positivereg","positiveregulator","possource","positivesource"]
#
#         #Get column header names as list
#         col_names = list(reading_df.columns)
#
#         #formate to bare minimum information
#         bare_cols = [x.lower().replace(" ","").replace("_","").replace("-","") for x in col_names]
#
#         #Check intersection of accepted column names and file column names
#         if {len(set(s_name_list) & set(bare_cols)) == 1 &
#             len(set(s_type_list) & set(bare_cols)) == 1 & len(set(s_id_list) & set(bare_cols)) == 1}:
#             #If minimum necessary columns are found, define variables for the column header
#             pos_source_name = col_names[bare_cols.index((set(s_name_list) & set(bare_cols)).pop())]
#             pos_source_type = col_names[bare_cols.index((set(s_type_list) & set(bare_cols)).pop())]
#             pos_source_id = col_names[bare_cols.index((set(s_id_list) & set(bare_cols)).pop())]
#             neg_source_name = pos_source_name.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
#             neg_source_type = pos_source_type.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
#             neg_source_id = pos_source_id.replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
#             #store column header names in a dictionary
#             reading_cols = {"pos_source_name" : pos_source_name, "pos_source_type" : pos_source_type, "pos_source_id" : pos_source_id,
#                             "neg_source_name" : neg_source_name, "neg_source_type" : neg_source_type, "neg_source_id" : neg_source_id}
#             #Now for the attributes:
#             for x in atts:
#                 s_att_list = [pre + x.lower().replace(" ","") for pre in s_att_pre]
#                 #made sure the attribute is in the reading columns
#                 if len(set(s_att_list) & set(bare_cols)) == 1:
#                     #if it's found, add it to the reading columns
#                     #add the attribute for the target
#                     reading_cols['pos_source_'+x.lower().replace(" ","_")] = col_names[bare_cols.index((set(s_att_list) & set(bare_cols)).pop())]
#                     reading_cols['neg_source_'+x.lower().replace(" ","_")] = reading_cols['pos_source_'+x.lower().replace(" ","_")].replace("positive","negative").replace("pos","neg").replace("Positive","Negative").replace("Pos","Neg")
#                 else:
#                     raise ValueError("Attribute \""+x+"\" was not found in your LEE input document."+"\n"+
#                     "Please check your file and try again")
#         else:
#             raise ValueError("Your LEE input is missing information."+"\n"+
#             "VIOLIN requires the following information: Name, Type, and ID of target node and regulators")
#         #End relative column name retrieval
#
#         reading_df['Reg Sign'] = pd.Series().astype(object)
#         #Move everything to "Positive" columns, add regulation sign
#         for y in range(len(list(reading_df[reading_cols['pos_source_name']]))):
#             if reading_df.loc[y,reading_cols['pos_source_name']] == 'nan':
#                 for each in [z for z in list(reading_cols.keys()) if 'pos' in z]:
#                     reading_df.loc[y,reading_cols[each]] = reading_df.loc[y,reading_cols[each.replace('pos','neg')]]
#                     reading_df.loc[y,'Reg Sign'] = 'decreases'
#             else: reading_df.loc[y,'Reg Sign'] = 'increases'
#
#         #Delete "Negative" Columns
#         for every in [a for a in list(reading_cols.keys()) if 'neg' in a]:
#             reading_df = reading_df.drop(columns=[reading_cols[every]])
#         #Delete "Positive " from header columns
#         reading_df.columns = [s.replace('Positive ', '').replace('Pos','') for s in reading_df.columns]
#
#         return reading_df
#
#     else: raise ValueError("Unsupported action. This function takes /'separate/' or /'combine/' as action input")
