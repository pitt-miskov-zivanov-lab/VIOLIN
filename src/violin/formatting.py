"""
formatting.py

Handles the model and reading formatting functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
Updated May 2024 - Haomiao Luo
"""
import logging
import re
import httplib2 as http
import json
import time
import pandas as pd

pd.options.mode.chained_assignment = None

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

headers = {
    "Accept": 'application/json'
}

h = http.Http()

# define regex for valid characters in variable names
_VALID_CHARS = r'a-zA-Z0-9\_'

# valid element types
CANONICAL_TYPES = [
    'protein',
    'proteinfamily', 'proteincomplex', 'proteinfamilyproteincomplex',
    'rna', 'mrna', 'microrna', 'trna',  # RNA
    'chemical', 'simplechemical', 'chemicalcompound', 'compound',  # chemical
    'biologicalprocess', 'bioprocess',  # biological process
    # 'mutation', 'geneticmutation',  # mutation
]

CANONICAL_PROTEIN = [
    'protein', 'proteinfamily', 'proteincomplex', 'proteinfamilyproteincomplex'
]
CANONICAL_CHEMICAL = [
    'chemical', 'simplechemical', 'chemicalcompound', 'compound', 'chemicalfamily'
]
CANONICAL_BIOPROCESS = [
    'biologicalprocess', 'bioprocess'
]
CANONICAL_RNA = [
    'rna', 'mrna', 'microrna', 'trna'
]

_VAR_COL = 'Variable'
_IDX_COL = '#'
REQUIRED_MODEL = ['Element Name', 'Element Type', 'Element IDs', 'Variable', 'Positive Regulator List',
                  'Negative Regulator List']

TYPE_ABBR_DICT = {
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

SUBTYPE_ABBR_DICT = {
    "enzyme": "enz",
    "transporter": "tsp",
    "transcription-factor": "tsf",
    "transcription-repressor": "tsr",
    "transducer": "tsd",
    "kinase": "kin",
    "interferon": "ifn",
    "interleukin": "ilk",
    "subunit": "sub",
    "cytokine": "cyt",
    "tyrosine": "tyr",
    "receptor": "rec",
    "caspase": "cas",
    "phosphatase": "pho",
    "adaptor": 'ada',
    "peptidase": "pep",
    "cyclin": "cyc",
    "growth-factor": "gwf",
    "binding": "bin",
    "molecule": "mol",
    "oncogene": "onc",
    "proto-oncogene": "pnc",
    "suppressor": "sup",
    "tumor": "tum",
    "signaling": "sig",
    "biological": "bio",
    "process": "prc",
    "protein": "prt",
    "redox": "red",
    "metallopeptidase": "mtp",
    "nonhistone": "nhs",
    "nucleoprotein": "ncp",
    "hormone": "hor",
    "ligase": 'lgs',
    "ligand": "lgd",
    "regulator": "reg",
    "ubiquitin-protein": "ubp",
    "catalytic": "cat",
    "gtpase": "gtp",
    "reverse": "rvs",
    "transcriptase": "tst",
    "dehydrogenase": "dhy",
    "hydrogenase": "hyd",
    "peroxidase": "pox",
    "oxidase": "oxi",
    "glycoprotein": "glp",
    "necrosis-factor": "nec",
    "apoptosis": "apo",
    "active": "act"
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

    # Convert reading to lower case, to prevent issues with case difference
    reading = reading_df.apply(lambda x: x.astype(str).str.lower())
    # The columns that aren't used to determine duplicates (such as Paper ID or Evidence Text)
    remainder = [x for x in reading_df.columns if x not in col_names]

    # As VIOLIN Identifies duplicates, it merges attributes from the remainder list into a single cell
    # This is how we count the number of times an LEE appears, and keep track of paper IDs and evidence text
    counted_reading = reading.groupby(col_names)[remainder[0]].apply(list).reset_index(name=remainder[0])
    for x in range(1, len(remainder)):
        sub = reading.groupby(col_names)[remainder[x]].apply(list).reset_index(name=remainder[x])
        counted_reading[remainder[x]] = sub[remainder[x]]

    # Counting the number of duplicates
    counted_reading['Evidence Score'] = counted_reading[remainder[0]].str.len()
    # counted_reading.to_csv("Trying.csv")

    return counted_reading


def add_regulator_names_id(model_df):
    """
    This function converts the model regulator lists from BioRECIPE variables to the common element names and database
    identifiers

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
    # removes the initial values from the model dataframe, as they're not needed
    # Also adds new columns for the positive and negative regulator names and IDs
    col_headers = list(model_df.columns)
    model_df = model_df[col_headers]
    reg_col_list = ['Positive Regulator List', 'Negative Regulator List']
    model_df[reg_col_list] = model_df[reg_col_list].apply(lambda x: x.astype(str).str.lower())
    # Columns for positive
    model_df['Positive Names'] = pd.Series().astype(object)
    model_df['Positive IDs'] = pd.Series().astype(object)
    model_df['Negative Names'] = pd.Series().astype(object)
    model_df['Negative IDs'] = pd.Series().astype(object)

    # Convert Regulators
    for sign in ['Negative', 'Positive']:
        for y in range(model_df.shape[0]):
            if model_df[sign + ' Regulator List'][y] in ['', "nan"]:
                model_df.at[y, sign + ' Names'] = "nan"
                model_df.at[y, sign + ' IDs'] = "nan"

            else:
                reg_name = model_df[sign + " Regulator List"][y].split(",")
                if '' in reg_name:
                    reg_name.remove('')
                reg_id = []
                reg_var = reg_name.copy()
                model_df.at[y, sign + ' Regulator List'] = reg_var

                # find index for regulator in variable column, and copy the Element Name and IDs to the new columns
                for element in reg_name:
                    idx = list(model_df['Listname']).index(element)
                    reg_name[reg_name.index(element)] = model_df['Element Name'][idx]
                    # idx = list(model_df["Element Name"]).index(element)
                    # Since there are multiple IDs for each element, need to keep track of which
                    # IDs go with which regulator
                    reg_id.append(model_df["Element IDs"][idx])

                model_df.at[y, sign + ' Names'] = reg_name
                model_df.at[y, sign + ' IDs'] = reg_id

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

    # collect invalid element names in a list, so they can be removed everywhere in the model
    # find invalid characters in element names and names starting with numbers
    invalid_names = [
        x for x in model[_VAR_COL]
        if re.search(r'(^[0-9]+)', x.strip()) or re.search(r'([^' + _VALID_CHARS + ']+)', x.strip())
    ]

    if len(invalid_names) > 0:
        logging.info('Formatting variable names: ')

    # remove invalid characters at the start of the variable name
    replace_names = [re.sub(r'^[^' + _VALID_CHARS + ']+', '', x) for x in invalid_names]
    # replace invalid characters elsewhere in variable names
    replace_names = [re.sub(r'[^' + _VALID_CHARS + ']+', '_', x) for x in replace_names]

    # add ELE_ at the beginning of names starting with numbers
    replace_names = [re.sub(r'(^[0-9]+)', 'ELE_\\1', x) for x in replace_names]

    name_pairs = zip(invalid_names, replace_names)

    for (invalid_name, replace_name) in name_pairs:
        logging.info('%s -> %s' % (invalid_name, replace_name))
        model.replace(re.escape(invalid_name), re.escape(replace_name), regex=True, inplace=True)

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

    global CANONICAL_TYPES
    global CANONICAL_PROTEIN
    global CANONICAL_CHEMICAL
    global CANONICAL_BIOPROCESS
    global CANONICAL_RNA
    input_type = ''.join(re.findall(r'[A-z]+', input_type)) if str(input_type) != 'nan' else 'other'
    if input_type in CANONICAL_TYPES:
        if input_type in CANONICAL_PROTEIN:
            return 'protein'
        elif input_type in CANONICAL_CHEMICAL:
            return 'chemical'
        elif input_type in CANONICAL_BIOPROCESS:
            return 'bioprocess'
        elif input_type in CANONICAL_RNA:
            return 'rna'
    else:
        return input_type


def get_hgnc_symbol(hgnc_id, hgnc_dict=None, url='https://rest.genenames.org/fetch/hgnc_id'):
    if hgnc_dict==None:
        hgnc_dict = {}
    else: 
        pass
    # Check if we can find the number
    if '.' in hgnc_id:
        hgnc_id = hgnc_id.split('.')[0]
        pass
    elif bool(re.fullmatch(r'd+', hgnc_id)):
        pass
    else: 
        raise NotImplementedError(f"Cannot handle {hgnc_id}...")
    if hgnc_id in hgnc_dict.keys():
        symbol = hgnc_dict[hgnc_id]
    else:
        response, content = h.request(
            url + f'/{hgnc_id}',
            'GET',
            '',
            headers
        )
        data = json.loads(content)
        status_code = False; i = 0; symbol = ''
        while status_code != True and i < 10:
            try:
                i += 1
                response, content = h.request(
                    url + f'/{hgnc_id}',
                    'GET',
                    '',
                    headers)
    
                if response['status'] == '200':
                    symbol = data['response']['docs'][0]['symbol']
                    hgnc_dict[hgnc_id] = symbol
                    status_code = True
                else:
                    pass
            except Exception as e:
                print(e)
                time.sleep(1)

    return symbol


def get_element(reg_rule, layer):
    """
    This function parses the regulation rule and disentangle the symbol operators
    converting rule to a list of regulators

    Parameters
    ----------
    reg_rule: str, list
        A BioRECIPE Regulation Rule
    layer: int
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
                assert (layer == 0)
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
                        parentheses += 1
                    elif char == '}':
                        parentheses -= 1

                    if parentheses == 0:
                        cutpoint = index
                        break

                necessary_element = reg_element[1: cutpoint]
                enhence_element = reg_element[cutpoint + 2:-1]

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
                _list = [element for ele_list in split_comma_out_parentheses(reg_element[1:-1])
                         for element in get_element(ele_list, 1)]
                regulator_list += _list
            else:
                assert (',' not in reg_element)

                if reg_element[-1] == '^':
                    regulator_list.append(reg_element[0:-1])
                elif '&' in reg_element:
                    regulator_list.append(reg_element[1:-1])
                elif '*' in reg_element:
                    multiply_reg_list = reg_element.split('*')
                    for reg_ in multiply_reg_list:
                        if re.search(r'^[0-9]', reg_):
                            pass
                        elif not re.search(r'[a-zA-Z0-9_!]+', reg_):
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
            reg_list.append(reg_rule[start:index + 1])
        elif char == '(' or char == '{' or char == '[':
            parentheses += 1
        elif char == ')' or char == '}' or char == ']':
            parentheses -= 1
        elif char == ',' and parentheses == 0:
            reg_list.append(reg_rule[start:index])
            start = index + 1
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
    if str(model_df.loc[idx, 'Element Type']).replace(' ', '') not in TYPE_ABBR_DICT:
        ele_type = model_df.loc[idx, 'Element Type'].replace(' ', '')
    else:
        ele_type = TYPE_ABBR_DICT[model_df.loc[idx, 'Element Type'].replace(' ', '')]
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
    # FIXME: Only get first subtype (TBD for the other subtypes)
    subtype = subtype.split(',')[0]
    if subtype not in ['', 'nan']:
        for x in [subname for subname in subtype.replace('(', ' ').replace(')', ' ').split(' ') if
                  subname not in ['', ' ']]:
            if x.lower() not in SUBTYPE_ABBR_DICT.keys():
                list_.append(x.strip())
            else:
                list_.append(SUBTYPE_ABBR_DICT[x.lower().strip()])
        abbr = ''.join(list_)
    else:
        abbr = 'nan'
    return abbr

# TODO: implement with functionality and integrate with BioRECIPE

# TODO: implement the functionality and integrate with BioRECIPE
