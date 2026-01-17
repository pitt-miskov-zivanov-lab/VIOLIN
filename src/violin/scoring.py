"""
scoring.py

Handles the Match Score, Kind Score, and Epistemic Value functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd
import networkx as nx
from violin.numeric import get_attributes, find_element, compare
from violin.network import path_finding
from typing import List, Union
import logging 

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Default kind score dict - categories
# KIND_DICT = {"strong corroboration" : 2,
#                 "empty attribute" : 1,
#                 "indirect interaction" : 1,
#                 "path corroboration" : 1,
#                 "hanging extension" : 40,
#                 "full extension" : 40,
#                 "internal extension" : 40,
#                 "specification" : 30,
#                 "dir contradiction" : 10,
#                 "sign contradiction" : 10,
#                 "att contradiction" : 10,
#                 "dir mismatch" : 20,
#                 "path mismatch" : 20,
#                 "self-regulation" : 20}

KIND_DICT_A = {"strong corroboration" : 2,
                "empty attribute" : 1,
                "indirect interaction" : 3,
                "path corroboration" : 5,
                "specification" : 7,
                "hanging extension" : 40,
                "full extension" : 39,
                "internal extension" : 38,
                "dir contradiction" : 11,
                "sign contradiction" : 10,
                "att contradiction" : 9,
                "dir mismatch" : 20,
                "path mismatch" : 19,
                "self-regulation" : 18}

KIND_DICT_B = {"strong corroboration" : 2,
                "empty attribute" : 1,
                "indirect interaction" : 3,
                "path corroboration" : 5,
                "specification" : 7,
                "hanging extension" : 40,
                "full extension" : 39,
                "internal extension" : 38,
                "dir contradiction" : 11,
                "sign contradiction" : 10,
                "att contradiction" : 9,
                "dir mismatch" : 20,
                "path mismatch" : 19,
                "self-regulation" : 18,
                "flagged4" : 17,
                "flagged5" : 16}

MATCH_DICT = {"source present" : 1,
                "target present" : 100,
                "both present" : 10,
                "neither present" : 0.1}

# Default attributes list is empty
atts_list = []


def match_score(x:int, reading_df: pd.DataFrame, model_df: pd.DataFrame, match_values: dict=None) -> int:
    """
    This function calculates the Match Score for an interaction from the reading.

    Parameters
    ----------
    x : int
        A row index of the dataframe of Interaction set (IS) to be scored
    reading_df : pd.DataFrame
        The reading dataframe
    model_df : pd.DataFrame
        The model dataframe
    match_values : dict
        Dictionary assigning Match Score values
        Default values found in MATCH_DICT

    Returns
    -------
    match : int
        Match Score value
    """
    global MATCH_DICT
    regulated = False
    regulator = False

    reg_sign = reading_df.loc[x, 'Sign']

    if match_values is None:
        match_values = MATCH_DICT

    # Search for regulated from reading in model
    if (find_element("name",
                     reading_df.loc[x, 'Regulated Name'],
                     reading_df.loc[x, 'Regulated Type'],
                     model_df) != -1 or
        find_element("hgnc",
                     reading_df.loc[x, 'Regulated HGNC Symbol'],
                     reading_df.loc[x, 'Regulated Type'],
                     model_df) != -1 or
        find_element("id",
                     reading_df.loc[x, 'Regulated ID'],
                     reading_df.loc[x, 'Regulated Type'],
                     model_df,
                     reading_df.loc[x, 'Regulated Database']) != -1 ):
        regulated = True

    # Search for regulator from reading in model
    if (find_element("name",
                     reading_df.loc[x, 'Regulator Name'],
                     reading_df.loc[x, 'Regulator Type'],
                     model_df) != -1 or
        find_element("hgnc",
                     reading_df.loc[x, 'Regulator HGNC Symbol'],
                     reading_df.loc[x, 'Regulator Type'],
                     model_df) != -1 or
        find_element("id",
                     reading_df.loc[x, 'Regulator ID'],
                     reading_df.loc[x, 'Regulator Type'],
                     model_df,
                     reading_df.loc[x, 'Regulator Database']) != -1 ):
        regulator = True

    # Scoring definition
    # Only the regulator in the model
    if not regulated and regulator: match = match_values['source present']
    # Only the regulated in the model
    elif regulated and not regulator: match = match_values['target present']
    # Both regulator and regulated in the model
    elif regulated and regulator: match = match_values['both present']
    # Neither present in the model
    else: match = match_values['neither present']

    return match

def _kind_score_to_model_int_id(row_idx: int,
                           query_kind: int, 
                           kinds_list: List[Union[str, int]], 
                           model_t_indices: List[int], 
                           model_s_indices: List[int],
                           counter_list: List[str],
                           search_point: int=0,
                           count_path:bool=False,
                           hits: list=None) -> str:
    """Function to record the interaction ID of the model that is classified. The function recursively records the unique IDs
    Parameters
    ----------
    row_idx: int
        The row index of the interaction to be classified.
    query_kind : int
        The kind score to be classified.
    kinds_list : list
        The list of kind scores.
    model_t_indices : list
        The model target indices.   
    model_s_indices : list
        The model source indices.
    counter_list: list
        A list to record the interactions that are identified as corroborated or contradicted interaction in model.
    search_point : int
        The point to start searching in the list of kind scores.
    count_path: bool
        An indicator whether to count the path or not. Default is False, because the path should not be an interaction in model.
    
    Returns
    -------
    str
        The interaction ID of the model that is classified.
    """
    # TODO: Implement path index finding for UI
    if count_path:
        raise NotImplementedError("Path finding is not implemented yet.")
    
    # Confirm that queried kind is in the list of kinds, starting from the srearch point
    if query_kind not in kinds_list[search_point:]:
        logger.info("%s: The queried kind %s is either not in the list of kinds or counted already." % (row_idx, query_kind))
        # Exit function
        return
    
    # Find the index of the queried kind in the list of kinds
    current_kind_idx = kinds_list.index(query_kind, search_point)

    # Find the corresponding source and target indices in the model
    source_found = model_s_indices[current_kind_idx % len(model_s_indices)]
    target_found = model_t_indices[current_kind_idx // len(model_s_indices)]
    
    # Check if the interaction ID is already in the counter
    if '%s+%s' % (source_found, target_found) not in counter_list:
        counter_list.append('%s+%s' % (source_found, target_found))
        if hits is not None:
            hits.append(query_kind)
            hits.append(target_found)
            hits.append(source_found)
        return 
    
    # If the interaction ID is already in the counter, check next kind is same as queried kind
    else:
        _kind_score_to_model_int_id(row_idx,
                               query_kind, 
                               kinds_list, 
                               model_t_indices, 
                               model_s_indices, 
                               counter_list, 
                               search_point=current_kind_idx + 1,
                               hits=hits)
        
        
    
def kind_score(x: int,
               model_df: pd.DataFrame,
               reading_df: pd.DataFrame,
               graph: nx.DiGraph,
               counter: dict,
               kind_values: dict=None,
               attributes: list=None,
               classify_scheme: str='1',
               mi_cxn: str='d') -> int:
    """
    This function calculates the Kind Score for an interaction in the Interactions Set (iIS).
    The kind score will be used to represent the subcategories.
    For further details, please find out in: https://www.biorxiv.org/content/10.1101/2024.07.21.604448v1.

    Parameters
    ----------
    x : int
        The row index for an iIS.
    model_df: pd.DataFrame
        The model dataframe
    reading_df : pd.DataFrame
        The reading dataframe.
    graph : nx.DiGraph
        A directed graph of the model,used when function calls path_finding module.
    counter: dict
        A dictionary to record the interactions that are identified as corroborated or contradicted interaction in model.
        Default value is None.
    kind_values : dict
        Dictionary assigning Kind Score values.
        Default values found in KIND_DICT_A and KIND_DICT_B.
    attributes : list
        A list of attributes compared between the model and the machine reading output.
        Default is None.
    classify_scheme: str
        The scheme of the classification ('1', '2', and '3').
        Default is '1'.
    mi_cxn : str
        What connection type should be assigned to model interactions if not available.
        Accepted values are "d" (direct) or "i" (indirect).
        Deafult is "d".

    Returns
    -------
    kind : int
        Kind Score, score value.
    """
    # Initialize the parameters
    global MATCH_DICT

    assert (classify_scheme in ['1', '2', '3'])

    if kind_values is None:
        if classify_scheme in ['1', '2']:
            kind_values = KIND_DICT_A

        elif classify_scheme == '3':
            kind_values = KIND_DICT_B
    else:
        pass

    if attributes is None:
        match_score = []

    assert (mi_cxn in ['d', 'i'])

    ### Finding interaction in Interactions Set (iIS) attributes ###
    # Finding iIS regulator sign
    signs = ['Negative', 'Positive']
    if reading_df.loc[x, 'Sign'].lower() in ['activate', 'positive', 'increase']: reg_sign = 'Positive'
    else: reg_sign = 'Negative'
    signs.remove(reg_sign)
    opp_sign = signs[0]

    # Finding iIS Connection Type (if not in iIS input, default to indirect, 'i')
    if 'Connection Type' in reading_df.columns: iis_cxn_type = reading_df.loc[x, 'Connection Type']
    else: iis_cxn_type = 'i'

    # Add full location information, if user want to compare location of the element
    if 'Regulated Compartment' in attributes and 'Regulated Compartment ID' not in attributes:
        attributes.insert(attributes.index('Regulated Compartment') + 1, 'Regulated Compartment ID')
    elif 'Regulated Compartment ID' in attributes and 'Regulated Compartment' not in attributes:
        attributes.insert(attributes.index('Regulated Compartment ID'), 'Regulated Compartment')
    else:
        pass

    if 'Regulator Compartment' in attributes and 'Regulator Compartment ID' not in attributes:
        attributes.insert(attributes.index('Regulator Compartment') + 1, 'Regulator Compartment ID')
    elif 'Regulator Compartment ID' in attributes and 'Regulator Compartment' not in attributes:
        attributes.insert(attributes.index('Regulator Compartment ID'), 'Regulator Compartment')
    else:
        pass

    # Create list for attributes (i.e., location attributes, context attributes, influence attributes)
    reading_atts = attributes

    # Finding iIS other attributes
    if len(attributes) > 0:
        # Attributes for iIS index 'x'
        reading_atts = {att: reading_df.loc[x, att] for att in reading_atts}
    else:
        reading_atts = {}

    # Comparing to model
    source_hgnc = find_element("hgnc",
                               reading_df.loc[x, 'Regulator HGNC Symbol'],
                               reading_df.loc[x, 'Regulator Type'],
                               model_df)
    source_name = find_element("name",
                               reading_df.loc[x, 'Regulator Name'],
                               reading_df.loc[x, 'Regulator Type'],
                               model_df)
    source_id = find_element("id",
                             reading_df.loc[x, 'Regulator ID'],
                             reading_df.loc[x, 'Regulator Type'],
                             model_df,
                             reading_df.loc[x, 'Regulator Database'])

    target_hgnc = find_element("hgnc",
                               reading_df.loc[x, 'Regulated HGNC Symbol'],
                               reading_df.loc[x, 'Regulated Type'],
                               model_df)
    target_name = find_element("name",
                               reading_df.loc[x, 'Regulated Name'],
                               reading_df.loc[x, 'Regulated Type'],
                               model_df)
    target_id = find_element("id",
                             reading_df.loc[x, 'Regulated ID'],
                             reading_df.loc[x, 'Regulated Type'],
                             model_df,
                             reading_df.loc[x, 'Regulated Database'])

    # Both regulator (source) and regulated (target) node found in the model
    if (source_name != -1 or source_hgnc != -1 or source_id != -1) and \
            (target_name != -1 or target_hgnc != -1 or target_id != -1):
        # Find indices of regulator element (target) in model
        #FIXME: TBD for order
        # Privilege: HGNC > Name > ID
        if source_hgnc != -1: model_s_indices = source_hgnc
        elif source_name != -1: model_s_indices = source_name
        else: model_s_indices = source_id

        if target_hgnc != -1: model_t_indices = target_hgnc
        elif target_name != -1: model_t_indices = target_name
        else: model_t_indices = target_id

        kinds = []
        # print(f"s:{model_s_indices}, t:{model_t_indices}")
        # Loop over each instance of the target and source in the model (since the same element may exist multiple status
        for t_idx in model_t_indices:
            # Regulator list in model
            model_s_list = model_df.loc[t_idx, reg_sign+' Regulator List']  \
                if model_df.loc[t_idx, reg_sign+' Regulator List'] != 'nan' else 'nan'
            # Regulator list of opposite sign
            model_s_opp = model_df.loc[t_idx, opp_sign+' Regulator List'] \
                if model_df.loc[t_idx, opp_sign+' Regulator List'] != 'nan' else 'nan'

            for s_idx in model_s_indices:

                source_listname = model_df.loc[s_idx, 'Listname']

                target_listname = model_df.loc[t_idx, 'Listname']

                # MI with match direction, match sign
                if str(model_s_list) != "nan" and source_listname in model_s_list.split(','):
                    # Index of regulator name within regulator list
                    s_index = model_s_list.split(',').index(source_listname)
                    # Finding index MI regulator variable
                    model_s_variable = model_df.loc[t_idx,reg_sign+' Regulator List'].split(',')[s_index]
                    #model_s_element = list(model_df['Variable']).index(model_s_variable)

                    # Find MI connection type
                    if (reg_sign+' Connection Type List') in model_df.columns.values.tolist() and \
                            all(cxn_type.lower().strip() in ['i', 'd'] for cxn_type in
                                                        model_df.loc[t_idx, reg_sign+' Connection Type List'].split(',')):
                        # Connection type
                        mi_cxn_type = model_df.loc[t_idx,reg_sign+' Connection Type List'].split(",")[s_index]
                    else: mi_cxn_type = mi_cxn

                    # List of model attributes to compare to reading attributes
                    model_atts = get_attributes(t_idx, s_idx, reg_sign, model_df, attributes)


                    # If iIS ="I" and MI = "I" or iIS = "D" and MI = "D": check attributes
                    if (iis_cxn_type == "i" and mi_cxn_type == "i") or (iis_cxn_type == "d" and mi_cxn_type != "i"):

                        compare_atts = compare(model_atts, reading_atts)
                        # Strong Corroboration - perfect match
                        if compare_atts == 0:
                            kinds.append(kind_values['strong corroboration'])
                        # Weak corroboration - the iIS presents less information than the model interaction
                        elif compare_atts == 1:
                            kinds.append(kind_values['empty attribute'])
                        # Specification - the iIS presents new information
                        elif compare_atts == 2:
                            kinds.append(kind_values['specification'])
                        # Contradiction - the iIS presents information that disputes the model interaction
                        elif compare_atts == 3:
                            kinds.append(kind_values['att contradiction'])

                    # If iIS = "D" and MI = "I"
                    elif iis_cxn_type == "d" and mi_cxn_type == "i":
                        compare_atts = compare(model_atts, reading_atts)
                        # If attributes are non-contradictory: iIS is a specification
                        if compare_atts in [0,1,2]: kinds.append(kind_values['specification'])
                        # Else: iIS is a contradiction
                        elif compare_atts == 3: kinds.append(kind_values['att contradiction'])

                    # If iIS ="I" and MI = "D":
                    elif iis_cxn_type == "i" and mi_cxn_type == "d":
                        compare_atts = compare(model_atts, reading_atts)
                        #If attributes are non-contradictory: iIS is a weak corroboration
                        if compare_atts in [0,1,2]: kinds.append(kind_values['indirect interaction'])
                        #Else: iIS is a contradiction
                        elif compare_atts == 3: kinds.append(kind_values['att contradiction'])

                # MI with Matched direction, Mismatched sign
                elif str(model_s_opp) != "nan" and source_listname in model_s_opp.split(','):
                    reg_index = model_df.loc[t_idx, opp_sign + " Regulator List"].split(',').index(
                        source_listname)
                    # Finding connection type
                    if (reg_sign+' Connection Type List') in model_df.columns.values.tolist() and \
                            all(cxn_type.lower().strip() in ['i', 'd'] for cxn_type in
                                                        model_df.loc[t_idx, opp_sign+' Connection Type List'].split(',')):
                        #Connection type
                        mi_cxn_type = model_df.loc[t_idx, opp_sign + ' Connection Type List'].split(",")[reg_index]
                    else: mi_cxn_type = mi_cxn
                    # If iIS = "I" and MI = "D"
                    if iis_cxn_type == "i" and mi_cxn_type != "i":
                        if classify_scheme in ['1', '2']:
                            kinds.append(kind_values['sign contradiction'])

                        elif classify_scheme == '3':
                            kinds.append(kind_values['flagged5'])

                    else:
                        #iIS is a Sign Contradiction, regardless of connection type
                        kinds.append(kind_values['sign contradiction'])

                # MI with Mismatched direction, Matched sign
                elif (model_df.loc[s_idx, reg_sign + " Regulator List"] != "nan" and target_listname
                      in model_df.loc[s_idx, reg_sign + " Regulator List"].split(',')):
                    reg_index = model_df.loc[s_idx, reg_sign + " Regulator List"].split(',').index(
                        target_listname)
                    # Finding connection type
                    if (reg_sign + ' Connection Type List') in model_df.columns.values.tolist() and \
                            all(cxn_type.lower().strip() in ['i', 'd'] for cxn_type in
                                                        model_df.loc[s_idx, reg_sign+' Connection Type List'].split(',')):
                        # Connection type
                        mi_cxn_type = model_df.loc[s_idx, reg_sign + ' Connection Type List'].split(",")[reg_index]
                    else:
                        mi_cxn_type = mi_cxn

                    # Finding index MI regulator variable
                    model_reg_variable = model_df.loc[s_idx, reg_sign + ' Regulator List'].split(',')[reg_index]
                    #model_reg_element = list(model_df['Variable']).index(model_reg_variable)
                    # List of model attributes to compare to reading attributes

                    model_atts = get_attributes(s_idx, t_idx, reg_sign, model_df, attributes)

                    # iIS = "I" and MI = "I"
                    if iis_cxn_type == "i" and mi_cxn_type == "i":
                        kinds.append(kind_values['dir contradiction'])

                    # iIS = "D" and MI = "D"
                    elif iis_cxn_type == "d" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        if classify_scheme in ['1', '2']:
                            # If the attributes are not contradictory - Flagged for manual review
                            if compare_atts in [0, 1, 2]:
                                kinds.append(kind_values['dir mismatch'])

                            # Else - Contradiction
                            else:
                                kinds.append(kind_values['dir contradiction'])
                        elif classify_scheme == '3':
                            kinds.append(kind_values['dir contradiction'])

                        else:
                            raise ValueError('Enter a right scheme number (1, 2, or 3).')

                    # iIS = "I" and MI = "D"
                    elif iis_cxn_type == "i" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        if classify_scheme in ['1', '2']:
                            # If the attributes are not contradictory - Flagged for manual review
                            if compare_atts in [0, 1, 2]:
                                kinds.append(kind_values['dir mismatch'])

                            # Else - Contradiction
                            else:
                                kinds.append(kind_values['dir contradiction'])
                        elif classify_scheme == '3':
                            if compare_atts in [0, 1, 2]:
                                kinds.append(kind_values['dir contradiction'])

                            else:
                                kinds.append(kind_values['flagged4'])

                        else:
                            raise ValueError('Enter a right scheme number (1, 2, or 3).')

                    # iIS = "D" and MI = "I"
                    elif iis_cxn_type == "d" and mi_cxn_type == "i":
                        kinds.append(kind_values['dir contradiction'])

                #MI with Mismatched direction, Mismatched sign
                elif (model_df.loc[s_idx,opp_sign+" Regulator List"] != "nan" and target_listname in
                      model_df.loc[s_idx,opp_sign+" Regulator List"].split(',')):
                    reg_index = model_df.loc[s_idx,opp_sign+" Regulator List"].split(',').index(target_listname)
                    #Finding connection type
                    if (opp_sign+' Connection Type List') in model_df.columns.values.tolist()and \
                            all(cxn_type.lower().strip() in ['i', 'd'] for cxn_type in
                                                        model_df.loc[s_idx, opp_sign+' Connection Type List'].split(',')):
                        mi_cxn_type = model_df.loc[s_idx,opp_sign+' Connection Type List'].split(",")[reg_index]
                    else: mi_cxn_type = mi_cxn

                    #Finding index MI regulator variable
                    model_reg_variable = model_df.loc[s_idx,opp_sign+' Regulator List'].split(',')[reg_index]
                    #model_reg_element = list(model_df['Variable']).index(model_reg_variable)

                    #List of model attributes to compare to reading attributes
                    model_atts = get_attributes(s_idx, t_idx, opp_sign, model_df, attributes)

                    # iIS = "D" and MI = "D"
                    if iis_cxn_type == "d" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        if classify_scheme in ['1', '2']:
                            #If the attributes are not contradictory - Flagged for manual review
                            if compare_atts in [0,1,2]:
                                kinds.append(kind_values['dir mismatch'])

                            #Else - Contradiction
                            else: kinds.append(kind_values['dir contradiction'])

                        elif classify_scheme == '3':
                            if compare_atts in [0,1,2]:
                                kinds.append(kind_values['dir contradiction'])

                            else:
                                kinds.append(kind_values['dir mismatch'])

                        else:
                            raise ValueError('Enter a right scheme (1, 2, or 3).')
                    # iIS = "D" and MI = "i"
                    elif iis_cxn_type == "d" and mi_cxn_type == "i": kinds.append(kind_values['dir contradiction'])
                    # iIS = "i" and MI = "D"
                    elif iis_cxn_type == "i" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        #If the attributes are not contradictory - Flagged for manual review
                        if compare_atts in [0, 1, 2]: kinds.append(kind_values['dir mismatch'])
                        #Else - Contradiction
                        else:
                            if classify_scheme in ['1', '2']:
                                kinds.append(kind_values['dir contradiction'])


                            elif classify_scheme == '3':
                                kinds.append(kind_values['flagged5'])

                    # iIS = "i" and MI = "i"
                    elif iis_cxn_type == "i" and mi_cxn_type == "i": kinds.append(kind_values['dir contradiction'])

                else:
                    # If there is a self-regulation (regulator is both target and source)
                    if t_idx == s_idx:
                        kind = kind_values['self-regulation']
                        hits = ["flagged", kind,
                            t_idx, s_idx]
                    # If model does not contain interaction - check for path
                    else:
                        kinds.append(path_finding(source_listname,target_listname,reg_sign,model_df,graph,kind_values,iis_cxn_type,reading_atts,attributes,classify_scheme))

        hits = []
        if len(kinds) == 1:
            kind = kinds[0]
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            # The tracker part is added after kind value is signed, to avoid corruption of the functionality of original version.
            else:
                if kind in [kind_values['strong corroboration'],
                            kind_values['empty attribute'],
                            kind_values['indirect interaction'],
                            kind_values['specification']]:
                    
                    _kind_score_to_model_int_id(
                        x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration'], hits=hits.append('corroboration')
                    )

                elif int(kind) in [kind_values['dir contradiction'],
                            kind_values['sign contradiction'],
                            kind_values['att contradiction']]:
                    if classify_scheme == '2':
                        # CS2 involves path for the category of contradiction, 
                        # For CS2, all the kind scores that are classified in `finding_path` function are converted to strings instead of int.
                        if type(kind) == str:
                            kind = int(kind)
                        else:
                            _kind_score_to_model_int_id(
                                x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction'], hits=hits.append('contradiction')
                            )

                    else:
                        _kind_score_to_model_int_id(
                            x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction'], hits=hits.append('contradiction')
                        )

                else:
                    pass

        # Strong Corroboration
        elif kind_values['strong corroboration'] in kinds:
            kind = kind_values['strong corroboration']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                _kind_score_to_model_int_id(
                    x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration'], hits=hits.append('corroboration')
                )

        # Weak Corroboration
        elif kind_values['empty attribute'] in kinds:
            kind = kind_values['empty attribute']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:

                _kind_score_to_model_int_id(
                    x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration'], hits=hits.append('corroboration')
                )

        elif kind_values['indirect interaction'] in kinds:
            kind = kind_values['indirect interaction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                _kind_score_to_model_int_id(
                    x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration'], hits=hits.append('corroboration')
                )

        elif kind_values['path corroboration'] in kinds:
            kind = kind_values['path corroboration']
        elif kind_values['specification'] in kinds:
            kind = kind_values['specification']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                _kind_score_to_model_int_id(
                    x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration'], hits=hits.append('corroboration')
                )

        # Contradiction
        elif kind_values['dir contradiction'] in kinds or str(kind_values['dir contradiction']) in kinds:
            kind = kind_values['dir contradiction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                if classify_scheme == '2':
                    kind_value = [x for x in kinds if
                                   x in [kind_values['dir contradiction'], str(kind_values['dir contradiction'])]]
                    for _ in kind_value:
                        if type(_) == str:
                            pass
                        else:
                            _kind_score_to_model_int_id(
                                x, _, kinds, model_t_indices, model_s_indices, counter['contradiction'], hits=hits.append('contradiction')
                            )

                            break
                else:
                    _kind_score_to_model_int_id(
                        x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction'], hits=hits.append('contradiction')
                    )

        elif kind_values['sign contradiction'] in kinds or str(kind_values['sign contradiction']) in kinds:
            kind = kind_values['sign contradiction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                if classify_scheme == '2':
                    kind_value = [x for x in kinds if
                                   x in [kind_values['sign contradiction'], str(kind_values['sign contradiction'])]]
                    for _ in kind_value:
                        if type(_) == str:
                            pass
                        else:
                            _kind_score_to_model_int_id(
                                x, _, kinds, model_t_indices, model_s_indices, counter['contradiction'], hits=hits.append('contradiction')
                            )
                            break
                else:
                    _kind_score_to_model_int_id(
                        x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction'], hits=hits.append('contradiction')
                    )

        elif kind_values['att contradiction'] in kinds or str(kind_values['att contradiction']) in kinds:
            kind = kind_values['att contradiction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                if classify_scheme == '2':
                    kind_value = [x for x in kinds if x in [kind_values['att contradiction'], str(kind_values['att contradiction'])]]
                    for _ in kind_value:
                        if type(_) == str:
                            pass
                        else:
                            _kind_score_to_model_int_id(
                                x, _, kinds, model_t_indices, model_s_indices, counter['contradiction'], hits=hits.append('contradiction')
                            )
                            break
                else:
                    _kind_score_to_model_int_id(
                        x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction'], hits=hits.append('contradiction')
                    )

        # Extensions
        elif kind_values['hanging extension'] in kinds:
            kind = kind_values['hanging extension']
        elif kind_values['internal extension'] in kinds:
            kind = kind_values['internal extension']
        elif kind_values['full extension'] in kinds:
            kind = kind_values['full extension']

        # Flagged
        elif kind_values['dir mismatch'] in kinds:
            kind = kind_values['dir mismatch']
        elif kind_values['path mismatch'] in kinds:
            kind = kind_values['path mismatch']
        elif kind_values['self-regulation'] in kinds:
            kind = kind_values['self-regulation']
        # check if the classify scheme is version 3
        elif 'flagged4' in kind_values:
            if kind_values['flagged4'] in kinds:
                kind = kind_values['flagged4']
        elif 'flagged5' in kind_values:
            if kind_values['flagged5'] in kinds:
                kind = kind_values['flagged5']
        else: pass

    # Both Extension - Both nodes from reading not in model
    elif (source_id == -1 and source_name == -1 and source_hgnc == -1) and (target_id == -1 and target_name == -1 and target_hgnc == -1):
        kind = kind_values['full extension']
    # Hanging Extension - One from reading not in model
    else: kind = kind_values['hanging extension']

    try:
        return kind, hits
    except UnboundLocalError: 
        print("unable to find hit...")
        print(kind)
        print(kinds)
        return kind, None


def epistemic_value(x: int,reading_df: pd.DataFrame) -> float:
    """
    Finds the epistemic value of the interaction in Interaction Set (IS) (when available).

    Parameters
    ----------
    x : int
        The row index for an iIS.
    reading_df : pd.DataFrame
        An IS dataframe.

    Returns
    -------
    e_value : float
        The Epistemic Value; if there is no Epistemic Value available for the reading, default is 1 for all interactions in IS.
    """

    if 'Epistemic Value' in reading_df.columns:
        e_value =  reading_df.loc[x,'Epistemic Value']
    else: e_value = 1

    return e_value


def score_reading(reading_df: pd.DataFrame,
                model_df: pd.DataFrame,
                graph:nx.DiGraph,
                counter: dict=None,
                kind_values: dict=None,
                match_values: dict=None,
                attributes: list=atts_list,
                classify_scheme: str='1',
                mi_cxn: str='d',
                include_hits=False) -> pd.DataFrame:
    """
    This function creates new columns for the Match Score, Kind Score, Epistemic Value, and Total Score.
    it calls scoring functions and stores the values in the approriate column.

    Parameters
    ----------
    reading_df : pd.DataFrame
        The reading dataframe.
    model_df : pd.DataFrame
        The model dataframe.
    graph : nx.DiGraph
        directed graph of the model, necessary for calling kind_score module.
    counter: dict
        A dictionary for counting the corrobrated and contradicted interaction.
        defulat value is None and ignore the counting step.
    kind_values : dict
        Dictionary assigning Kind Score values.
        Default values found in KIND_DICT_A and KIND_DICT_B.
    match_values : dict
        Dictionary assigning Match Score values.
        Default values found in MATCH_DICT.
    attributes : list
        List of attributes compared between the model and the machine reading output.
        Default is None.
    classify_scheme: str
        The scheme of the classification.
        Default value is '1'.
    include_hits: bool
        Whether to include the hits information

    Returns
    -------
    scored = reading_df : pd.DataFrame
        reading dataframe with added scores.
    """
    assert (classify_scheme in ['1', '2', '3'])
    assert (mi_cxn in ['d', 'i'])

    if kind_values is None:
        if classify_scheme in ['1', '2']:
            kind_values = KIND_DICT_A

        elif classify_scheme == '3':
            kind_values = KIND_DICT_B
    else:
        pass

    if match_values is None:
        match_values = MATCH_DICT


    #Create new DF columns for score calculations
    scored_reading_df = reading_df.copy()
    scored_reading_df['Match Score'] = pd.Series()
    scored_reading_df['Kind Score'] = pd.Series()
    scored_reading_df['Epistemic Value'] = pd.Series()
    scored_reading_df['Total Score'] = pd.Series()

    #Calculate scores
    collections = []
    for x in range(reading_df.shape[0]):
        scored_reading_df.at[x,'Match Score'] = match_score(x,reading_df,model_df, match_values)
        scored_reading_df.at[x,'Kind Score'], interaction_hit = kind_score(x,model_df,reading_df,graph, counter,kind_values,attributes,classify_scheme,mi_cxn)
        scored_reading_df.at[x,'Epistemic Value'] = epistemic_value(x,reading_df)
        scored_reading_df.at[x,'Total Score'] =  ((scored_reading_df.at[x,'Evidence Score']*scored_reading_df.at[x,'Match Score'])+int(scored_reading_df.at[x,'Kind Score']))*scored_reading_df.at[x,'Epistemic Value']
        if interaction_hit != None:
            collections.append({"category": interaction_hit[0], 
                            "kind_score": str(interaction_hit[1]),
                            "target_index": interaction_hit[2],
                            "source_index": interaction_hit[3],
                            "interaction_index": x})
    if include_hits:
        return scored_reading_df, collections
    else: 
        return scored_reading_df
