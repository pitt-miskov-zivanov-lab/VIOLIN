"""
numeric.py

Handles the element finding and comparison functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd

def get_attributes(A_idx, B_idx, sign, model_df, attrs: list, path=False):
    """
    The function get the attributes of the interaction in model
    Parameters
    ----------
    A_idx: int
    A represents the element
    B_idx: int
    B represents the regulator
    sign: string
        'positive' or 'negative'
    model_df: pd.DataFrame object
        DataFrame of BioRECIPE model file
    attrs: list
        attributes list for reading file
    path: boolean
        influence attributes will be empty if only path is found in model

    Returns
    -------
    model_atts, dictionary
    """
    model_attrs = {attr: x for attr, x in zip(attrs, ['nan'] * len(attrs))}
    # Check if user input redundant attributes
    if set(attrs).issubset({'Regulator Compartment', 'Regulator Compartment ID',
                            'Regulated Compartment', 'Regulated Compartment ID',
                            'Mechanism', 'Site',
                            'Cell Line', 'Cell Type', 'Tissue Type', 'Organism',
                            }):
        pass
    else:
        raise ValueError('VIOLIN accepts attributes '
                         'Regulator Compartment, Regulator Compartment ID,'
                         'Regulated Compartment, Regulated Compartment ID,'
                         'Mechanism, Site,'
                         'Cell Line, Cell Type, Tissue Type, Organism')

    # For influence attributes
    if path:  # influence attributes will be empty if only path is found in model
        pass
    else:
        assert (sign in ['Positive', 'Negative'])
        source_position = model_df.loc[A_idx, f'{sign} Regulator List'].split(',').index(
            model_df.loc[B_idx, 'Listname'])
        for a in ['Mechanism', 'Site']:
            if a in attrs:
                if model_df.at[A_idx, f'{sign} {a} List'] != 'nan' and \
                        model_df.loc[A_idx, f'{sign} {a} List'].split(',')[source_position] not in ['none',
                                                                                                    'nan', '']:
                    model_attrs[a] = model_df.loc[A_idx, f'{sign} {a} List'].split(',')[source_position]
                else:
                    pass
            else:
                pass

    # For context attributes
    for a in ['Cell Line', 'Cell Type', 'Tissue Type', 'Organism']:
        if a in attrs:
            if model_df.at[A_idx, a] not in ['none', 'nan', '']:
                model_attrs[a] = model_df.at[A_idx, a]
            else:
                pass
        else:
            pass

    # For element attributes
    if 'Regulated Compartment' in attrs:
        A_location = model_df.loc[A_idx, 'Compartment']
        A_location_id = model_df.loc[A_idx, 'Compartment ID']
        model_attrs['Regulated Compartment'] = A_location if A_location.lower() not in ['none', 'nan', ''] else 'nan'
        model_attrs['Regulated Compartment ID'] = A_location_id if A_location_id.lower() not in ['none', 'nan',
                                                                                                 ''] else 'nan'

    if 'Regulator Compartment' in attrs:
        B_location = model_df.loc[B_idx, 'Compartment']
        B_location_id = model_df.loc[B_idx, 'Compartment ID']
        model_attrs['Regulator Compartment'] = B_location if B_location.lower() not in ['none', 'nan', ''] else 'nan'
        model_attrs['Regulator Compartment ID'] = B_location_id if B_location_id.lower() not in ['none', 'nan',
                                                                                                 ''] else 'nan'

    return model_attrs


def find_element(search_type,
                 element_name,
                 element_type,
                 model_df,
                 id_db=None ):
    """
    This function finds the correct indices of an element within the model.
    Because elements can exist as multiple types (Protein, RNA, gene, etc.),
    this function checks the element name/ID along with the element type.
    Function may return a list, if a given element of a specific type 
    exists with varying attributes (such as different locations)

    Parameters
    ----------
    search_type : str
        Whether the element is being searched for by 'name' or 'ID'
    element_name: str
        The name (or ID) of the element being searched for
    id_db: str
        provide id's database
    element_type: str
        The type of element searched for ('protein', 'protein family', etc.)
    model_df: pd.DataFrame
        The model dataframe

    Returns
    -------
    location : list
        All rows of the model spreadsheet in which the element is found (returns -1 if not found)
    """

    # Searching for element by HGNC symbol
    if search_type == "hgnc":
        # indices of all instances of an element in the model
        indices = [i for i, x in enumerate(list(model_df['Element HGNC Symbol'])) if
                   element_name in x and element_name != 'nan']
    # Searching for element by name
    elif search_type == "name":
        # indices of all instances of an element in the model
        indices = [i for i, x in enumerate(list(model_df['Element Name'])) if
                   element_name in x and element_name != 'nan']
    # Searching for element by ID
    elif search_type == "id":
        # indices of all instances of an element in the model
        indices = [i for i, x in enumerate(list(model_df['Element IDs'])) if (element_name in x) \
                   and (element_name != 'nan') and (id_db == model_df.loc[i, 'Element Database'])]

    else:
        indices = []
    indices_list = []

    # Searching for matching element type
    for idx in indices:
        # if model_df.loc[idx, "Element Type"] == element_type:
        if model_df.loc[idx, "Element Type"] in element_type or element_type in model_df.loc[idx, "Element Type"]:
            indices_list.append(idx)

    # If element has been found, return a list of its locations within the model
    if len(indices_list) > 0:
        return indices_list
    # Value -1: means element not found
    else:
        return -1


def compare(model_atts: dict, reading_atts: dict):
    """
    Compares a list of model attributes to the corresponding LEE attributes, returns numeric value\n
        Attributes are the same (strong corroboration): 0 \n
        Some or all LEE attributes are missing (weak corroboration): 1\n
        Some or all of the model attributes are missing (specification): 2\n
        One or more model attribute differs from the LEE attributes (contradiction): 3

    Parameters
    ----------
    model_atts : dict
        dictionary of attributes from the model interaction
    reading_atts : dict
        dictionary of attributes from the literature extracted event (LEE)

    Returns
    -------
    value : int
        Numerical representation of comparison outcome
    """

    # First, individually compare each attribute between the machine reading output and model
    # Add outcome to list
    compare_atts = []

    s_location_atts, t_location_atts = [], []
    for model, reading in zip(model_atts.keys(), reading_atts.keys()):
        # Check compartment information
        if reading in ['Regulated Compartment', 'Regulated Compartment ID']:

            if model_atts[model] == reading_atts[reading]:
                t_location_atts.append(0)
            elif (model_atts[model] == "nan") and (reading_atts[reading] == "nan"):
                t_location_atts.append(0)
            # Reading attribute is not available
            elif (model_atts[model] != "nan") and (reading_atts[reading] == "nan"):
                t_location_atts.append(1)
            # Model attribute is not available
            elif (model_atts[model] == "nan") and (reading_atts[reading] != "nan"):
                t_location_atts.append(2)
            # Else: both model and reading attributes are available, but they differ
            else:
                t_location_atts.append(3)

        elif reading in ['Regulator Compartment', 'Regulator Compartment ID']:

            if model_atts[model] == reading_atts[reading]:
                s_location_atts.append(0)
            elif (model_atts[model] == "nan") and (reading_atts[reading] == "nan"):
                s_location_atts.append(0)
            # Reading attribute is not available
            elif (model_atts[model] != "nan") and (reading_atts[reading] == "nan"):
                s_location_atts.append(1)
            # Model attribute is not available
            elif (model_atts[model] == "nan") and (reading_atts[reading] != "nan"):
                s_location_atts.append(2)
            # Else: both model and reading attributes are available, but they differ
            else:
                s_location_atts.append(3)

        # Check other attributes
        else:
            if model_atts[model] == reading_atts[reading]:
                compare_atts.append(0)
            elif (model_atts[model] == "nan") and (reading_atts[reading] == "nan"):
                compare_atts.append(0)
            # Reading attribute is not available
            elif (model_atts[model] != "nan") and (reading_atts[reading] == "nan"):
                compare_atts.append(1)
            # Model attribute is not available
            elif (model_atts[model] == "nan") and (reading_atts[reading] != "nan"):
                compare_atts.append(2)
            # Else: both model and reading attributes are available, but they differ
            else:
                compare_atts.append(3)

    s_location_atts = [0] * len(s_location_atts) if any(x == 0 for x in s_location_atts) else s_location_atts
    t_location_atts = [0] * len(t_location_atts) if any(x == 0 for x in t_location_atts) else t_location_atts

    compare_atts += s_location_atts
    compare_atts += t_location_atts

    # Final outcome determined by the set of comparisons in compare_atts[]
    # Strong Corroboration - perfect match over all attributes
    if 0 in compare_atts and len(set(compare_atts)) == 1:
        value = 0

    # Weak Corroboration - model contains more information than the reading for all attributes
    elif 1 in compare_atts and len(set(compare_atts)) == 1:
        value = 1
    # Weak Corroboration - attributes either match or the model contains more information than the LEE
    elif 0 in compare_atts and 1 in compare_atts and len(set(compare_atts)) == 2:
        value = 1

    # Specification - LEE contains one or more attributes the model does not
    elif 2 in compare_atts and len(set(compare_atts)) == 1:
        value = 2
    # Specification - attributes either match or the LEE contains additional attributes
    elif 0 in compare_atts and 2 in compare_atts and len(set(compare_atts)) == 2:
        value = 2
    # Specification - model or LEE contains more information, depending on attribute
    elif 1 in compare_atts and 2 in compare_atts and len(set(compare_atts)) == 2:
        value = 2

    # Contradiction - some or all of the attributes differ
    elif 3 in compare_atts:
        value = 3

    # Specification - combination of perfect match, model contains more information, reading contains more information
    else:
        value = 2

    return value
