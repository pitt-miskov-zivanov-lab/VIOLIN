import pandas as pd
import networkx as nx
from typing import List, Union
from violin.formatting import get_type
from violin.numeric import find_element, get_attributes, compare
from violin.network import path_finding
from violin.utils import ElementMatcher
import logging

element_matcher = ElementMatcher()

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

_LABEL2SUB = {
    'strong corroboration': 'corroboration',
    'empty attribute': 'corroboration',
    'indirect interaction': 'corroboration',
    'path corroboration': 'corroboration',
    'specification': 'corroboration',
    'full extension': 'extension',
    'hanging extension': 'extension',
    'internal extension': 'extension',
    'dir contradiction': 'contradiction',
    'sign contradiction': 'contradiction',
    'att contradiction': 'contradiction',
    'dir mismatch': 'flagged',
    'path mismatch': 'flagged',
    'self-regulation': 'flagged',
    'flagged4': 'flagged',
    'flagged5': 'flagged',
}

MATCH_DICT = {"source present" : 1,
                "target present" : 100,
                "both present" : 10,
                "neither present" : 0.1}

_SCORE2LABEL_A = {v: k for k, v in KIND_DICT_A.items()}
_SCORE2LABEL_B = {v: k for k, v in KIND_DICT_B.items()}

# Default attributes list is empty
atts_list = []


def _kind_score_to_model_int_id(row_idx: int,
                           query_kind: int, 
                           kinds_list: List[Union[str, int]], 
                           model_t_indices: List[int], 
                           model_s_indices: List[int],
                           counter_list: List[str],
                           search_point: int=0,
                           count_path:bool=False) -> str:
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
        return 
    
    # If the interaction ID is already in the counter, check next kind is same as queried kind
    else:
        _kind_score_to_model_int_id(row_idx,
                               query_kind, 
                               kinds_list, 
                               model_t_indices, 
                               model_s_indices, 
                               counter_list, 
                               search_point=current_kind_idx + 1)

def kind_score(x: int,
               model_df: pd.DataFrame,
               reading_df: pd.DataFrame,
               graph: nx.DiGraph,
               counter: dict,
               kind_values: dict=None,
               attributes: list=None,
               classify_scheme: str='1',
               mi_cxn: str='d',
               name_match_threshold: float=None,
               match_sim_metric: str='jaccard') -> int:
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
    name_match_threshold : float
        The threshold of the match between the reading entity and the model entity
        e.g., 0.8 means that at least the similarity of compared reading entity string and model entity string
        reach 0.8, the entity pair will be considered as a match. 
    match_sim_metric : str
        The similarity metric for calculating the similarity.
        Available options: "jaccard", "edit_sim", and "grounding". Default is "jaccard".

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
            _SCORE2LABEL = _SCORE2LABEL_A

        elif classify_scheme == '3':
            kind_values = KIND_DICT_B
            _SCORE2LABEL = _SCORE2LABEL_B
    else:
        _SCORE2LABEL = {v: k for k, v in kind_values.items()}
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

    if match_sim_metric == 'grounding':
        if element_matcher.model_df is None:
            element_matcher.initialize_model_curies(model_df)
        source_grd, source_grd_confs =  element_matcher.find_element(
            reading_df.loc[x, 'Regulator Name'], 
            threshold=name_match_threshold,
        )
        target_grd, target_grd_confs = element_matcher.find_element(
            reading_df.loc[x, 'Regulated Name'], 
            threshold=name_match_threshold,
        )
        
    else:
        source_hgnc, source_hgnc_confs = find_element("hgnc",
                                reading_df.loc[x, 'Regulator HGNC Symbol'],
                                reading_df.loc[x, 'Regulator Type'],
                                model_df,
                                threshold=name_match_threshold,
                                metric=match_sim_metric)
        source_name, source_name_confs = find_element("name",
                                reading_df.loc[x, 'Regulator Name'],
                                reading_df.loc[x, 'Regulator Type'],
                                model_df,
                                threshold=name_match_threshold,
                                metric=match_sim_metric)
        source_id, source_id_confs = find_element("id",
                                reading_df.loc[x, 'Regulator ID'],
                                reading_df.loc[x, 'Regulator Type'],
                                model_df,
                                reading_df.loc[x, 'Regulator Database'],
                                threshold=name_match_threshold,
                                metric=match_sim_metric)
        

        target_hgnc, target_hgnc_confs = find_element("hgnc",
                                reading_df.loc[x, 'Regulated HGNC Symbol'],
                                reading_df.loc[x, 'Regulated Type'],
                                model_df,
                                threshold=name_match_threshold,
                                metric=match_sim_metric)
        target_name, target_name_confs = find_element("name",
                                reading_df.loc[x, 'Regulated Name'],
                                reading_df.loc[x, 'Regulated Type'],
                                model_df,
                                threshold=name_match_threshold,
                                metric=match_sim_metric)
        target_id, target_id_confs = find_element("id",
                                reading_df.loc[x, 'Regulated ID'],
                                reading_df.loc[x, 'Regulated Type'],
                                model_df,
                                reading_df.loc[x, 'Regulated Database'],
                                threshold=name_match_threshold,
                                metric=match_sim_metric)

    # Both regulator (source) and regulated (target) node found in the model
    hits = []
    if match_sim_metric == 'grounding':
        has_source = source_grd != -1
        has_target = target_grd != -1
    else:
        has_source = any(x != -1 for x in [source_hgnc, source_name, source_id])
        has_target = any(x != -1 for x in [target_hgnc, target_name, target_id])
    
    
    if has_source and has_target:
        # Find indices of regulator element (target) in model
        if match_sim_metric == 'grounding':
            model_s_indices = source_grd
            model_s_confs = source_grd_confs
            model_t_indices = target_grd
            model_t_confs = target_grd_confs
        else:
            #FIXME: TBD for order
            # Privilege: HGNC > Name > ID
            if source_hgnc != -1: 
                model_s_indices = source_hgnc
                model_s_confs = source_hgnc_confs
            elif source_name != -1: 
                model_s_indices = source_name
                model_s_confs = source_name_confs
            else: 
                model_s_indices = source_id
                model_s_confs = source_id_confs

            if target_hgnc != -1: 
                model_t_indices = target_hgnc
                model_t_confs = target_hgnc_confs
            elif target_name != -1: 
                model_t_indices = target_name
                model_t_confs = target_name_confs
            else: 
                model_t_indices = target_id
                model_t_confs = target_id_confs

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

                source_variable = model_df.loc[s_idx, 'Variable']

                target_variable = model_df.loc[t_idx, 'Variable']

                # MI with match direction, match sign
                if str(model_s_list) != "nan" and source_variable in model_s_list.split(','):
                    # Index of regulator name within regulator list
                    s_index = model_s_list.split(',').index(source_variable)
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
                elif str(model_s_opp) != "nan" and source_variable in model_s_opp.split(','):
                    reg_index = model_df.loc[t_idx, opp_sign + " Regulator List"].split(',').index(
                        source_variable)
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
                elif (model_df.loc[s_idx, reg_sign + " Regulator List"] != "nan" and target_variable
                      in model_df.loc[s_idx, reg_sign + " Regulator List"].split(',')):
                    reg_index = model_df.loc[s_idx, reg_sign + " Regulator List"].split(',').index(
                        target_variable)
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
                elif (model_df.loc[s_idx,opp_sign+" Regulator List"] != "nan" and target_variable in
                      model_df.loc[s_idx,opp_sign+" Regulator List"].split(',')):
                    reg_index = model_df.loc[s_idx,opp_sign+" Regulator List"].split(',').index(target_variable)
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
                        t_list_idx= model_t_indices.index(t_idx)
                        s_list_idx = model_s_indices.index(s_idx)

                        hits = ["flagged", kind,
                            t_idx, s_idx, model_t_confs[t_list_idx], model_s_confs[s_list_idx]]
                    # If model does not contain interaction - check for path
                    else:
                        kinds.append(path_finding(source_variable,target_variable,reg_sign,model_df,graph,kind_values,iis_cxn_type,reading_atts,attributes,classify_scheme))
        
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
                        x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration']
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
                                x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction']
                            )

                    else:
                        _kind_score_to_model_int_id(
                            x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction']
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
                    x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration']
                )

        # Weak Corroboration
        elif kind_values['empty attribute'] in kinds:
            kind = kind_values['empty attribute']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:

                _kind_score_to_model_int_id(
                    x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration']
                )

        elif kind_values['indirect interaction'] in kinds:
            kind = kind_values['indirect interaction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                _kind_score_to_model_int_id(
                    x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration']
                )

        elif kind_values['path corroboration'] in kinds:
            kind = kind_values['path corroboration']
            # Pseudo-hit for path corroboration

        elif kind_values['specification'] in kinds:
            kind = kind_values['specification']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                _kind_score_to_model_int_id(
                    x, kind, kinds, model_t_indices, model_s_indices, counter['corroboration']
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
                                x, _, kinds, model_t_indices, model_s_indices, counter['contradiction']
                            )

                            break
                else:
                    _kind_score_to_model_int_id(
                        x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction']
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
                                x, _, kinds, model_t_indices, model_s_indices, counter['contradiction']
                            )
                            break
                else:
                    _kind_score_to_model_int_id(
                        x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction']
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
                                x, _, kinds, model_t_indices, model_s_indices, counter['contradiction']
                            )
                            break
                else:
                    _kind_score_to_model_int_id(
                        x, kind, kinds, model_t_indices, model_s_indices, counter['contradiction']
                    )

        # Extensions
        # TODO: this branch will never be taken, deprecated in the future version
        elif kind_values['hanging extension'] in kinds:
            kind = kind_values['hanging extension']
            if counter is None:
                pass
            else:
                raise NotImplementedError('hanging extension counter is not implemented yet.')
            
        elif kind_values['internal extension'] in kinds:
            kind = kind_values['internal extension']
        
        #TODO: this branch will never be taken, deprecated in the future version
        elif kind_values['full extension'] in kinds:
            kind = kind_values['full extension']
            if counter is None:
                pass
            else:
                raise NotImplementedError('full extension counter is not implemented yet.')

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
    elif not has_source and not has_target:
        kind = kind_values['full extension']
    # Hanging Extension - One from reading not in model
    else: 
        kind = kind_values['hanging extension']


    # Handle hits corresponding to kind
    # Since kind is assigned based on priority, 
    # the rest of the condition statements are regardless of the priority order.
    if kind == kind_values['hanging extension']:
        if has_source:
            assert has_target == False
            if match_sim_metric == 'grounding':
                source_idx = source_grd[0]
                regulator_conf = source_grd_confs[0]
            else:
                source_idx = source_hgnc[0] if source_hgnc != -1 else source_name[0] if source_name != -1 else source_id[0]
                regulator_conf = source_hgnc_confs[0] if source_hgnc != -1 else source_name_confs[0] if source_name != -1 else source_id_confs[0]
            hits = ['extension', kind, '', source_idx, 0.0, regulator_conf]

        if has_target:
            assert has_source == False
            if match_sim_metric == 'grounding':
                target_idx = target_grd[0]
                regulated_conf = target_grd_confs[0]
            else:
                target_idx = target_hgnc[0] if target_hgnc != -1 else target_name[0] if target_name != -1 else target_id[0]
                regulated_conf = target_hgnc_confs[0] if target_hgnc != -1 else target_name_confs[0] if target_name != -1 else target_id_confs[0]
            hits = ['extension', kind, target_idx, '', regulated_conf, 0.0]
    
    elif kind == kind_values['full extension']:
        hits = ['extension', kind, '', '', 0.0, 0.0]
    
    elif set(['flagged4', 'flagged5']).intersection(set(kind_values.keys())) and kind in [kind_values['flagged4'], kind_values['flagged5']]:
        raise NotImplementedError('Flagged4 and Flagged5 hits are not implemented yet.')
    
    elif kind == kind_values['self-regulation']:
        pass
    
    else:
        t_list_idx = kinds.index(kind) // len(model_s_indices)
        s_list_idx = kinds.index(kind) % len(model_s_indices)
        hits = [_LABEL2SUB[_SCORE2LABEL[kind]], kind, 
                model_t_indices[t_list_idx],
                model_s_indices[s_list_idx],
                model_t_confs[t_list_idx],
                model_s_confs[s_list_idx]]


    # invalid hits
    if len(hits) != 6:
        return kind, None
    return kind, hits


