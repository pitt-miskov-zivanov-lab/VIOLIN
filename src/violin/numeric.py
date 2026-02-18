"""
numeric.py

Handles the element finding and comparison functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
"""
import re
import pandas as pd
from typing import Tuple, Union, List
from violin.utils import get_similarity_score


def get_attributes(A_idx: int, B_idx: int, sign: str, model_df: pd.DataFrame, attrs: list, path: bool=False) -> dict:
    """
    The function gets the attributes of the interaction in model, available attributes includes
    [Regulator Compartment, Regulator Compartment ID, Regulated Compartment, Regulated Compartment ID,
    Mechansim, Site, Cell Line, Cell Type, Tissue Type, Organism]. If Regulator Compartment is selected,
    Regulator Compartment ID will also be selected.

    Parameters
    ----------
    A_idx: int
        A row index of element A in the input model dataframe.
    B_idx: int
        A row index of element B in the input model dataframe.
    sign: str
        A sign of the interaction, available options: 'positive' or 'negative'.
    model_df: pd.DataFrame
        A DataFrame of a model with BioRECIPE format.
    attrs: list
        An attributes list for interactions file.
    path: bool
        An indicator if it is path interaction. Attributes will be empty if only path is found in model.

    Returns
    -------
    model_atts: dict
        An dict of attributes for a model interaction.
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
        raise ValueError('acceptable attributes '
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
            model_df.loc[B_idx, 'Variable'])
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


def find_element(search_type: str,
                 element_name: str,
                 element_type: str,
                 model_df: pd.DataFrame,
                 id_db: str=None,
                 threshold: float=None,
                 metric: str="jaccard") -> Union[Tuple[List[int], List[float]], Tuple[int, list]]:
    """
    This function finds the correct indices of an element within the model.
    Because elements can exist as multiple types (protein, RNA, gene, etc.),
    this function checks the element name/ID along with the element type.
    Function may return a list, if a given element of a specific type
    exists with varying attributes (such as different locations).

    Parameters
    ----------
    search_type : str
        An identifier of the element, available options are 'hgnc', 'name', and 'id'.
    element_name: str
        A name (or ID) of the element being searched for.
    element_type: str
        A type of element ('protein', 'protein family', etc.)
    model_df: pd.DataFrame
        A model dataframe within BioRECIPE format.
    id_db: str
        A database name for provided identifier.
    threshold: float
        A similarity score threshold for element matching, only used for 'name', default is None,
        which means the losest match will be accepted.
    metric: str
        A metric for calculating the similarity score, available options are 'jaccard' and 'edit_sim', 'None'.
        If 'None', the similarity score will not be calculated

    Returns
    -------
    location : list|int
        All row indices of the model spreadsheet in which the element is found (returns -1 if not found).
    """

    if metric == "none": 
        # Threshold should not be provided when metric is None.
        threshold = None

    # Searching for element by HGNC symbol
    if search_type == "hgnc":
        # indices of all instances of an element in the model
        sim_score = model_df['Element HGNC Symbol'].apply(get_similarity_score, name2=element_name, metric=metric)

        # Default option
        if not threshold:
            indices = [i for i, x in enumerate(list(model_df['Element HGNC Symbol'])) if
                    element_name in x and element_name != 'nan']
        else:
            indices = sim_score[sim_score >= threshold].index.tolist()

        conf_scores = [sim_score[i] for i in indices]
    # Searching for element by name
    elif search_type == "name":

        sim_score = model_df['Element Name'].apply(get_similarity_score, name2=element_name, metric=metric)

        # Default option
        if not threshold:
            # indices of all instances of an element in the model
            indices = [i for i, x in enumerate(list(model_df['Element Name'])) if
                    element_name in x and element_name != 'nan']
        else: 
            # take into account similarity score
            indices = sim_score[sim_score >= threshold].index.tolist()
        conf_scores = [sim_score[i] for i in indices]
    # Searching for element by ID, calculating the similarity score of between a ID and IDs list is not meaningful,
    # so we only check if the provided ID is a substring of the model element IDs,
    # and the model element database matches the provided database name. 
    elif search_type == "id":
        # TODO: support for more accurate ID confidence score calculation
        elements_ids = model_df['Element IDs'].apply(
            lambda x: re.split(r'[\/\|\,\;\ ]+', x) if x != 'nan' else []
        )

        # [model_len, ids_list_len]
        all_sim_scores = elements_ids.apply(lambda ids_list: [get_similarity_score(element_name, model_id, metric=metric) \
                                                   for model_id in ids_list] if ids_list else [])

        # TODO: remove the legacy code for ID matching
        # indices of all instances of an element in the model
        # indices = [i for i, x in enumerate(list(model_df['Element IDs'])) if (element_name in x) \
        #         and (element_name != 'nan') and (id_db == model_df.loc[i, 'Element Database'])]
        # TODO end

        # find exactly match IDs list 
        # [match_len, ]
        indices = [i for i, x in enumerate(elements_ids) if (element_name in x) \
                and (element_name != 'nan') and (id_db == model_df.loc[i, 'Element Database'])]
        # [match_len, ]
        match_position = [elements_ids[i].index(element_name) for i in indices]
        # take out all similarity score of the matched IDs
        # [match_len, ]
        match_sim_indices = {i: pos for i, pos in zip(indices, match_position)}
        # [match_len, ]
        match_sim_scores = [all_sim_scores[i][pos] for i, pos in match_sim_indices.items()]
        if not threshold:
            pass
        else:
            # Remove indices with similarity score below the threshold
            # [filtered_indices_len, ]
            indices = [x for i, x in enumerate(indices) if match_sim_scores[i] >= threshold]
            
        # Take out confidence scores of the remaining indices
        # [filtered_indices_len, ]
        conf_scores = [all_sim_scores[i][match_sim_indices[i]] for i in indices]

    else:
        indices = []
    indices_list = []

    # Searching for matching element type
    for idx in indices:
        if model_df.loc[idx, "Element Type"] == element_type:
            indices_list.append(idx)

    # If element has been found, return a list of its locations within the model
    if len(indices_list) > 0:
        return indices_list, conf_scores
    # Value -1: means element not found
    else:
        return -1, []


def compare(model_atts: dict, reading_atts: dict) -> int:
    """
    Compares a list of model attributes to the corresponding interaction attributes, returns numeric value
        - Attributes are the same (strong corroboration): 0
        - Some or all LEE attributes are missing (weak corroboration): 1
        - Some or all of the model attributes are missing (specification): 2
        - One or more model attribute differs from the LEE attributes (contradiction): 3

    Parameters
    ----------
    model_atts : dict
        A dictionary of attributes for a model interaction.
    reading_atts : dict
        A dictionary of attributes for an event from a literactions list.

    Returns
    -------
    value : int
        The numerical representation of comparison outcome.
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
    # Check if compartment ID is matched or not
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


