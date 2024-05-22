import pandas as pd
import re
import logging
import datetime
from collections import OrderedDict
# define regex for valid characters in variable names
_VALID_CHARS = r'a-zA-Z0-9\_'

# valid element types
_VALID_TYPES = [
    'protein', 'protein family', 'protein complex',
    'rna', 'mrna', 'gene', 'chemical', 'biological process'
    ]

_VAR_COL = 'Variable'
_IDX_COL = '#'

BIORECIPE_MODEL = OrderedDict([("number_row", "#"),
                                    ("element_name", "Element Name"),
                                    ("element_type", "Element Type"),
                                    ("element_subtype", "Element Subtype"),
                                    ("element_hgnc_symbol", "Element HGNC Symbol"),
                                    ("element_database", "Element Database"),
                                    ("element_ids", "Element IDs"),
                                    ("compartment", "Compartment"),
                                    ("compartment_id", "Compartment ID"),
                                    ("cell_line", "Cell Line"),
                                    ("cell_type", "Cell Type"),
                                    ("tissue_type", "Tissue Type"),
                                    ("organism", "Organism"),
                                    ("positive_regulator_list", "Positive Regulator List"),
                                    ("positive_connection_type_list", "Positive Connection Type List"),
                                    ("positive_mechanism_list", "Positive Mechanism List"),
                                    ("positive_site_list", "Positive Site List"),
                                    ("negative_regulator_list", "Negative Regulator List"),
                                    ("negative_connection_type_list", "Negative Connection Type List"),
                                    ("negative_mechanism_list", "Negative Mechanism List"),
                                    ("negative_site_list", "Negative Site List"),
                                    ("score_list", "Score List"),
                                    ("statements_list", "Statements List"),
                                    ("paper_ids_list", "Paper IDs List"),
                                    ("positive_regulation_rule", "Positive Regulation Rule"),
                                    ("negative_regulation_rule", "Negative Regulation Rule"),
                                    ("variable", "Variable"),
                                    ("value_type","Value Type"),
                                    ("levels", "Levels"),
                                    ("state_list0", "State List 0"),
                                    ("state_list1", "State List 1"),
                                    ("state_list2", "State List 2"),
                                    ("const_off", "Const OFF"),
                                    ("const_on", "Const ON"),
                                    ("increment", "Increment"),
                                    ("spontaneous", "Spontaneous"),
                                    ("balancing", "Balancing"),
                                    ("delay", "Delay"),
                                    ("update_group", "Update Group"),
                                    ("update_rate", "Update Rate"),
                                    ("update_rank", "Update Rank")
                                   ])


VIOLIN_MODEL = OrderedDict([("number_row", "#"),
                                    ("element_name", "Element Name"),
                                    ("element_type", "Element Type"),
                                    ("element_subtype", "Element Sub-type"),
                                    ("element_ids", "Element IDs"),
                                    ("cell_line", "Cell Line"),
                                    ("cell_type", "Cell Type"),
                                    ("organism", "Organism"),
                                    ("tissue_type", "Tissue Type"),
                                    ("compartment", "Location"),
                                    ("compartment_id", "Location ID"),
                                    ("element_note", "Element NOTES"),
                                    ("blank_column1", "1"),
                                    ("variable", "Variable"),
                                    ("blank_column2", "2"),
                                    ("blank_column3", "3"),
                                    ("positive_regulator_list", "Positive Regulators"),
                                    ("positive_connection_type_list", "Positive Connection Type"),
                                    ("negative_regulator_list", "Negative Regulators"),
                                    ("negative_connection_type_list", "Negative Connection Type")
                                  ])


def get_model(model_file):

    """
    This function loads model into a DataFrame and standardize column names

    Parameters
    ----------
    model_file: str
        The name of the model file

    Returns
    -------
    model: pd.DataFrame
        A dataframe containing the model
    """

    global _VALID_CHARS
    global _VAR_COL
    global _IDX_COL

    index_col_name = _IDX_COL
    var_col_name = _VAR_COL
    pos_reg_col_name = 'Positive Regulation Rule'
    pos_list_col_name = 'Positive List'
    neg_reg_col_name = 'Negative Regulation Rule'
    neg_list_col_name = 'Negative List'
    reg_list_col_name = 'Regulators'
    element_name_col_name = 'Element Name'
    ids_col_name = 'Element IDs'
    type_col_name = 'Element Type'

    # Load the input file containing elements and regulators
    model_sheets = pd.ExcelFile(model_file)
    # get the model from the first sheet, will check the other sheets for truth tables later
    model = model_sheets.parse(0,na_values='NaN',keep_default_na=False,index_col=None)

    # check model format
    if 'element attributes' in [x.lower() for x in model.columns]:
        # drop two header rows and set column names to third row
        model = model.rename(columns=model.iloc[1]).drop([0,1]).set_index(index_col_name)

    # get other sheets
    # TODO: parse truth tables here? or just return other sheets separately?
    if len(model_sheets.sheet_names) > 1:
        df_other_sheets = {sheet : model_sheets.parse(sheet,na_values='NaN',keep_default_na=False) \
            for sheet in model_sheets.sheet_names[1:]}
    else:
        df_other_sheets = ''

    # format model columns
    input_col_X = [
            x.strip() for x in model.columns
            if ('variable' in x.lower())
            ]
    input_col_A = [
            x.strip() for x in model.columns
            if ('positive regulation rule' in x.lower())
            ]
    input_col_I = [
            x.strip() for x in model.columns
            if ('negative regulation rule' in x.lower())
            ]
    input_col_initial = [
            x.strip() for x in model.columns
            if ('state list' in x.lower())
            ]

    input_col_name = [
            x.strip() for x in model.columns
            if ('element name' in x.lower())
            ]
    input_col_ids = [
            x.strip() for x in model.columns
            if ('element ids' in x.lower())
            ]
    input_col_type = [
            x.strip() for x in model.columns
            if ('element type' in x.lower())
            ]

    # check for all required columns or duplicate colummns
    if (len(input_col_X) == 0
            or len(input_col_A) == 0
            or len(input_col_I) == 0
            or len(input_col_initial) == 0
            ):
        raise ValueError(
                'Missing one or more required columns in input file: '
                'Variable, Positive Regulation Rule, Negative Regulation Rule, State List'
                )
    elif (len(input_col_X) > 1
            or len(input_col_A) > 1
            or len(input_col_I) > 1
            ):
        raise ValueError('Duplicate column of: Variable, Positive Regulation Rule, Negative Regulation Rule')

    if (len(input_col_name) == 0
            or len(input_col_ids) == 0
            or len(input_col_type) == 0
            ):
        raise ValueError(
                'Missing one or more required column names: '
                'Element Name, Element IDs, Element Type'
                )
    elif (len(input_col_name) > 1
            or len(input_col_ids) > 1
            or len(input_col_type) > 1
            ):
        raise ValueError(
                'Duplicate column of: '
                'Element Name, Element IDs, Element Type'
                )

    # TODO: check for other columns here as they are needed

    # processing
    # use # column or index to preserve order of elements in the model
    if index_col_name in model.columns:
        model.set_index(index_col_name,inplace=True)


    # model = model.reset_index()
    # standardize column names
    model = model.rename(
        index=str,
        columns={
            'index': index_col_name,
            input_col_X[0]: var_col_name,
            input_col_A[0]: pos_reg_col_name,
            input_col_I[0]: neg_reg_col_name,
            input_col_name[0]: element_name_col_name,
            input_col_ids[0]: ids_col_name,
            input_col_type[0]: type_col_name
        })

    # format invalid variable names
    model = format_variable_names(model)

    # standardize element types
    model['Element Type'] = model['Element Type'].apply(get_type)

    # set variable name as the index
    model.set_index(var_col_name,inplace=True)

    # check for empty indices
    if '' in model.index:
        raise ValueError('Missing variable names')
        # model = model.drop([''])

    # parse regulation functions into lists of regulators
    model[pos_list_col_name] = model[pos_reg_col_name].apply(
            lambda x: [y.strip() for y in re.findall('['+_VALID_CHARS+']+',x)]
            )
    model[neg_list_col_name] = model[neg_reg_col_name].apply(
            lambda x: [y.strip() for y in re.findall('['+_VALID_CHARS+']+',x)]
            )
    model[reg_list_col_name] = model.apply(
            lambda x:
            set(list(x[pos_list_col_name]) + list(x[neg_list_col_name])),
            axis=1
            )

    model.fillna('',inplace=True)

    return model


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


def norm_model(model_file, save_dir):

    """
    This function converts a BioRECIPE model to VIOLIN, regulator lists are got from regulation rules
        # model_file: model filename
        # save_dir: the directory of VIOLIN format model

    Parameters
    ----------
    model_file: str
        model filename
    save_dir: str
        directory of VIOLIN format model
    Returns
    -------

    """

    # Read file as BioRECIPE format
    model_df = get_model(model_file)
    # check if the regulator and regulation columns are empty or not
    c = 0
    for sign in ['Positive', 'Negative']:
        regulator = all(x in ['', 'Nan', 'nan'] for x in model_df[f'{sign} Regulator List'])
        regulation = all(x in ['', 'Nan', 'nan'] for x in model_df[f'{sign} Regulation Rule'])

        if regulator and regulation:
            c+=1
    if c > 1:
        raise ValueError(
            "The regulation rule and list columns are both empty, please fill at least one column out")
    model_df = model_df.fillna('')

    # Register Model columns
    VIOLIN_cols = list(VIOLIN_MODEL.values())
    VIOLIN_attr = list(VIOLIN_MODEL.keys())
    output_df = pd.DataFrame(columns=VIOLIN_cols)
    variable_index = model_df.index

    for row in range(len(model_df)):
        input_index = variable_index[row]
        for sign in ['Positive', 'Negative']:
            # if regulator and c <= 1:
            #     if model_df.loc[input_index, f'{sign} Regulation Rule'] == '':
            #         model_df.loc[input_index, f'{sign} Regulator List'] = ''
            #     else:
            #         model_df.loc[input_index, f'{sign} Regulator List'] = ','.join(list(get_element(model_df.loc[input_index, f'{sign} Regulation Rule'], 0)))
            # elif not regulator:
            #     reg_list = model_df.loc[input_index, f'{sign} Regulator List'].split(',')
            #     model_df.loc[input_index, f'{sign} Regulator List'] =','.join(list(reg_list))
            #
            # else: pass
            if model_df.loc[input_index, f'{sign} Regulation Rule'] == '':
                model_df.loc[input_index, f'{sign} Regulator List'] = ''
            else:
                model_df.loc[input_index, f'{sign} Regulator List'] = ','.join(list(get_element(model_df.loc[input_index, f'{sign} Regulation Rule'], 0)))

        for key, value in BIORECIPE_MODEL.items():
            if value not in model_df:
                pass
            else:
                if key == "number_row":
                    output_df.loc[input_index, '#'] = row if model_df.loc[input_index, "Element Type"] not in ['', 'nan', 'NaN'] else 'X'
                elif key in VIOLIN_attr:
                    output_df.loc[input_index, VIOLIN_MODEL[key]] = model_df.loc[input_index, BIORECIPE_MODEL[key]]
                else:
                    pass
        output_df.loc[input_index, 'Variable'] = input_index

    df = pd.DataFrame(output_df, columns=VIOLIN_cols)
    df.to_excel(save_dir, index=False)
    return


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
