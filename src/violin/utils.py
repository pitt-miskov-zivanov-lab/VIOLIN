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
