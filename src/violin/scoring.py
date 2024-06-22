"""
scoring.py

Handles the Match Score, Kind Score, and Epistemic Value functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd
from src.violin.numeric import get_attributes, find_element, compare
from src.violin.network import path_finding
from src.violin.formatting import get_listname

kind_dict = {"strong corroboration" : 2, 
                "empty attribute" : 1,
                "indirect interaction" : 1,
                "path corroboration" : 1,
                "hanging extension" : 40, 
                "full extension" : 40, 
                "internal extension" : 40, 
                "specification" : 30, 
                "dir contradiction" : 10,
                "sign contradiction" : 10,
                "att contradiction" : 10,
                "dir mismatch" : 20,
                "path mismatch" : 20,
                "self-regulation" : 20}
match_dict = {"source present" : 1, 
                "target present" : 100, 
                "both present" : 10, 
                "neither present" : 0.1}

# Default attributes list is empty
atts_list = []

def match_score(x, reading_df, model_df, embedding_match, match_values = match_dict):
    """
    This function calculates the Match Score for an interaction from the reading

    Parameters
    ----------
    x : int
        The line of the reading dataframe with the interaction to be scored
    reading_df : pd.DataFrame
        The reading dataframe
    model_df : pd.DataFrame
        The model dataframe
    match_values : dict
        Dictionary assigning Match Score values
        Default values found in match_dict

    Returns
    -------
    match : int
        Match Score value
    """
    regulated = False
    regulator = False
    reg_sign = reading_df.loc[x, 'Sign']

    # Search for regulated from reading in model
    if (find_element("name",
                     reading_df.loc[x, 'Regulated Name'],
                     reading_df.loc[x, 'Regulated Type'],
                     model_df,
                     embedding_match) != -1 or
        find_element("hgnc",
                     reading_df.loc[x, 'Regulated HGNC Symbol'],
                     reading_df.loc[x, 'Regulated Type'],
                     model_df,
                     embedding_match) != -1 or
        find_element("id",
                     reading_df.loc[x, 'Regulated ID'],
                     reading_df.loc[x, 'Regulated Type'],
                     model_df,
                     embedding_match) != -1 ):
        regulated = True

    # Search for regulator from reading in model
    if (find_element("name",
                     reading_df.loc[x, 'Regulator Name'],
                     reading_df.loc[x, 'Regulator Type'],
                     model_df,
                     embedding_match) != -1 or
        find_element("hgnc",
                     reading_df.loc[x, 'Regulator HGNC Symbol'],
                     reading_df.loc[x, 'Regulator Type'],
                     model_df,
                     embedding_match) != -1 or
        find_element("id",
                     reading_df.loc[x, 'Regulator ID'],
                     reading_df.loc[x, 'Regulator Type'],
                     model_df,
                     embedding_match) != -1 ):
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


def kind_score(x,
               model_df,
               reading_df,
               graph,
               embedding_match,
               counter,
               kind_values = kind_dict,
               attributes = atts_list,
               classify_scheme = '1',
               mi_cxn = 'd'):
    """
    This function calculates the Kind Score for an interaction in the reading

    Parameters
    ----------
    x : int
        The line of the reading dataframe with the interaction to be scored
    model_df: pd.DataFrame
        The model dataframe
    reading_df : pd.DataFrame
        The reading dataframe
    graph : nx.DiGraph
        directed graph of the model,used when function calls path_finding module
    counter: dict
        A dictionary to record the interactions that are identified as corroborated or contradicted interaction in model
        default value is None
    kind_values : dict
        Dictionary assigning Kind Score values
        Default values found in kind_dict
    attributes : list
        List of attributes compared between the model and the machine reading output
        Default is None
    classify_scheme: str
        The scheme of the classification ('1', '2', and '3')
        Default is '1'
    mi_cxn : str
        What connection type should be assigned to model interactions if not available
        Accepted values are "d" (direct) or "i" (indirect)
        Deafult is "d"

    Returns
    -------
    kind : int
        Kind Score score value
    """

    ### Finding LEE attributes ###
    # Finding LEE regulator sign
    signs = ['Negative', 'Positive']
    if reading_df.loc[x, 'Sign'].lower() in ['activate', 'positive', 'increase']: reg_sign = 'Positive'
    else: reg_sign = 'Negative'
    signs.remove(reg_sign)
    opp_sign = signs[0]

    # Finding LEE Connection Type (if not in lee input, default to indirect, 'i')
    if 'Connection Type' in reading_df.columns: lee_cxn_type = reading_df.loc[x, 'Connection Type']
    else: lee_cxn_type = 'i'

    # Add full location information, if user want to compare location of the element
    if 'Regulated Compartment' in attributes and 'Regulated Compartment ID' not in attributes:
        attributes.insert(attributes.index('Regulated Compartment') + 1, 'Regulated Compartment ID')
    elif 'Regulated Compartment ID' in attributes and 'Regulated Compartment' not in attributes:
        attributes.insert(attributes.index('Regulated Compartment ID'), 'Regulated Compartment')

    if 'Regulator Compartment' in attributes and 'Regulator Compartment ID' not in attributes:
        attributes.insert(attributes.index('Regulator Compartment') + 1, 'Regulator Compartment ID')
    elif 'Regulator Compartment ID' in attributes and 'Regulator Compartment' not in attributes:
        attributes.insert(attributes.index('Regulator Compartment ID'), 'Regulator Compartment')
    else:
        pass

    # Create list for attributes (i.e., location attributes, context attributes, influence attributes)
    reading_atts = attributes

    # Finding Lee other attributes
    if len(attributes) > 0:
        # Attributes for LEE index 'x'
        reading_atts = {att: reading_df.loc[x, att] for att in reading_atts}
    else:
        reading_atts = {}

    # Comparing to model
    source_name = find_element("name",
                               reading_df.loc[x, 'Regulator Name'],
                               reading_df.loc[x, 'Regulator Type'],
                               model_df,
                               embedding_match)
    source_hgnc = find_element("hgnc",
                               reading_df.loc[x, 'Regulator HGNC Symbol'],
                               reading_df.loc[x, 'Regulator Type'],
                               model_df,
                               embedding_match)
    source_id = find_element("id",
                             reading_df.loc[x, 'Regulator ID'],
                             reading_df.loc[x, 'Regulator Type'],
                             model_df,
                             embedding_match)

    target_name = find_element("name",
                               reading_df.loc[x, 'Regulated Name'],
                               reading_df.loc[x, 'Regulated Type'],
                               model_df,
                               embedding_match)
    target_hgnc = find_element("hgnc",
                               reading_df.loc[x, 'Regulated HGNC Symbol'],
                               reading_df.loc[x, 'Regulated Type'],
                               model_df,
                               embedding_match)
    target_id = find_element("id",
                             reading_df.loc[x, 'Regulated ID'],
                             reading_df.loc[x, 'Regulated Type'],
                             model_df,
                             embedding_match)

    # Both regulator (source) and regulated (target) node found in the model
    if (source_name != -1 or source_hgnc != -1 or source_id != -1) and \
            (target_name != -1 or target_hgnc != -1 or target_id != -1):
        # Find indices of regulator element (target) in model
        #FIXME: TBD for order
        # Priority: HGNC > Name > ID
        if source_hgnc != -1: model_s_indices = source_hgnc
        elif source_name != -1: model_s_indices = source_name
        else: model_s_indices = source_id

        if target_hgnc != -1: model_t_indices = target_hgnc
        elif target_name != -1: model_t_indices = target_name
        else: model_t_indices = target_id

        kinds = []

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


                    # If LEE ="I" and MI = "I" or LEE = "D" and MI = "D": check attributes
                    if (lee_cxn_type == "i" and mi_cxn_type == "i") or (lee_cxn_type == "d" and mi_cxn_type != "i"):

                        compare_atts = compare(model_atts, reading_atts)
                        # Strong Corroboration - perfect match
                        if compare_atts == 0:
                            kinds.append(kind_values['strong corroboration'])
                        # Weak corroboration - the LEE presents less information than the model interaction
                        elif compare_atts == 1:
                            kinds.append(kind_values['empty attribute'])
                        # Specification - the LEE presents new information
                        elif compare_atts == 2:
                            kinds.append(kind_values['specification'])
                        # Contradiction - the LEE presents information that disputes the model interaction
                        elif compare_atts == 3:
                            kinds.append(kind_values['att contradiction'])

                    # If LEE = "D" and MI = "I"
                    elif lee_cxn_type == "d" and mi_cxn_type == "i":
                        compare_atts = compare(model_atts, reading_atts)
                        # If attributes are non-contradictory: LEE is a specification
                        if compare_atts in [0,1,2]: kinds.append(kind_values['specification'])
                        # Else: LEE is a contradiction
                        elif compare_atts == 3: kinds.append(kind_values['att contradiction'])

                    # If LEE ="I" and MI = "D":
                    elif lee_cxn_type == "i" and mi_cxn_type == "d":
                        compare_atts = compare(model_atts, reading_atts)
                        #If attributes are non-contradictory: LEE is a weak corroboration
                        if compare_atts in [0,1,2]: kinds.append(kind_values['indirect interaction'])
                        #Else: LEE is a contradiction
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
                    # If LEE = "I" and MI = "D"
                    if lee_cxn_type == "i" and mi_cxn_type != "i":
                        if classify_scheme in ['1', '2']:
                            kinds.append(kind_values['sign contradiction'])
                        elif classify_scheme == '3':
                            kinds.append(kind_values['flagged5'])
                    else:
                        #LEE is a Sign Contradiction, regardless of connection type
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

                    # LEE = "I" and MI = "I"
                    if lee_cxn_type == "i" and mi_cxn_type == "i":
                        kinds.append(kind_values['dir contradiction'])

                    # LEE = "D" and MI = "D"
                    elif lee_cxn_type == "d" and mi_cxn_type != "i":
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

                    # LEE = "I" and MI = "D"
                    elif lee_cxn_type == "i" and mi_cxn_type != "i":
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

                    # LEE = "D" and MI = "I"
                    elif lee_cxn_type == "d" and mi_cxn_type == "i":
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

                    # LEE = "D" and MI = "D"
                    if lee_cxn_type == "d" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        if classify_scheme in ['1', '2']:
                            #If the attributes are not contradictory - Flagged for manual review
                            if compare_atts in [0,1,2]: kinds.append(kind_values['dir mismatch'])
                            #Else - Contradiction
                            else: kinds.append(kind_values['dir contradiction'])

                        elif classify_scheme == '3':
                            if compare_atts in [0,1,2]: kinds.append(kind_values['dir contradiction'])
                            else: kinds.append(kind_values['dir mismatch'])
                        else:
                            raise ValueError('Enter a right scheme (1, 2, or 3).')
                    # LEE = "D" and MI = "i"
                    elif lee_cxn_type == "d" and mi_cxn_type == "i": kinds.append(kind_values['dir contradiction'])
                    # LEE = "i" and MI = "D"
                    elif lee_cxn_type == "i" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        #If the attributes are not contradictory - Flagged for manual review
                        if compare_atts in [0, 1, 2]: kinds.append(kind_values['dir mismatch'])
                        #Else - Contradiction
                        else:
                            if classify_scheme in ['1', '2']: kinds.append(kind_values['dir contradiction'])

                            elif classify_scheme == '3':
                                kinds.append(kind_values['flagged5'])
                    # LEE = "D" and MI = "D"
                    elif lee_cxn_type == "i" and mi_cxn_type == "i": kinds.append(kind_values['dir contradiction'])

                else:
                    # If there is a self-regulation (regulator is both target and source)
                    if t_idx == s_idx:
                        kind = kind_values['self-regulation']
                    # If model does not contain interaction - check for path
                    else:
                        kinds.append(path_finding(source_listname,target_listname,reg_sign,model_df,graph,kind_values,lee_cxn_type,reading_atts,attributes,classify_scheme))

        if len(kinds) == 1:
            kind = kinds[0]
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                if kind in [kind_values['strong corroboration'],
                            kind_values['empty attribute'],
                            kind_values['indirect interaction'],
                            kind_values['specification']]:
                    counter['corroboration'].append('{}+{}'.format(
                        model_t_indices[kinds.index(kind) // len(model_s_indices)],
                        model_s_indices[kinds.index(kind) % len(model_s_indices)])
                    )

                elif kind in [kind_values['dir contradiction'],
                            kind_values['sign contradiction'],
                            kind_values['att contradiction']]:
                    counter['contradiction'].append('{}+{}'.format(
                        model_t_indices[kinds.index(kind) // len(model_s_indices)],
                        model_s_indices[kinds.index(kind) % len(model_s_indices)])
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
                counter['corroboration'].append('{}+{}'.format(
                    model_t_indices[kinds.index(kind) // len(model_s_indices)],
                    model_s_indices[kinds.index(kind) % len(model_s_indices)])
                )

        # Weak Corroboration
        elif kind_values['empty attribute'] in kinds:
            kind = kind_values['empty attribute']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                counter['corroboration'].append('{}+{}'.format(
                    model_t_indices[kinds.index(kind) // len(model_s_indices)],
                    model_s_indices[kinds.index(kind) % len(model_s_indices)])
                )
        elif kind_values['indirect interaction'] in kinds:
            kind = kind_values['indirect interaction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                counter['corroboration'].append('{}+{}'.format(
                    model_t_indices[kinds.index(kind) // len(model_s_indices)],
                    model_s_indices[kinds.index(kind) % len(model_s_indices)])
                )
        elif kind_values['path corroboration'] in kinds:
            kind = kind_values['path corroboration']
        elif kind_values['specification'] in kinds:
            kind = kind_values['specification']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                counter['corroboration'].append('{}+{}'.format(
                    model_t_indices[kinds.index(kind) // len(model_s_indices)],
                    model_s_indices[kinds.index(kind) % len(model_s_indices)])
                )

        # Contradiction
        elif kind_values['dir contradiction'] in kinds:
            kind = kind_values['dir contradiction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                counter['contradiction'].append('{}+{}'.format(
                    model_t_indices[kinds.index(kind) // len(model_s_indices)],
                    model_s_indices[kinds.index(kind) % len(model_s_indices)])
                )
        elif kind_values['sign contradiction'] in kinds:
            kind = kind_values['sign contradiction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                counter['contradiction'].append('{}+{}'.format(
                    model_t_indices[kinds.index(kind) // len(model_s_indices)],
                    model_s_indices[kinds.index(kind) % len(model_s_indices)])
                )
        elif kind_values['att contradiction'] in kinds:
            kind = kind_values['att contradiction']
            if counter is None:
                pass
            # Track every matched interaction that is classified as corroborated interaction or contradicted interaction
            else:
                counter['contradiction'].append('{}+{}'.format(
                    model_t_indices[kinds.index(kind) // len(model_s_indices)],
                    model_s_indices[kinds.index(kind) % len(model_s_indices)])
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

    return kind


def epistemic_value(x,reading_df):
    """
    Finds the epistemic value of the LEE (when available)

    Parameters
    ----------
    x : int
        The line of the reading dataframe with the interaction to be scored
    reading_df : pd.DataFrame
        The reading dataframe

    Returns
    -------
    e_value : float
        The Epistemic Value; if there is no Epistemic Value available for the reading, default is 1 for all LEEs
    """

    if 'Epistemic Value' in reading_df.columns:
        e_value =  reading_df.loc[x,'Epistemic Value']
    else: e_value = 1

    return e_value


def score_reading(reading_df, model_df, graph,
                  embedding_match=False, counter=None,
                  kind_values = kind_dict, match_values = match_dict,
                  attributes = atts_list, classify_scheme = '1', mi_cxn = 'd'):
    """
    Creates new columns for the Match Score, Kind Score, Epistemic Value, and Total Score.
    Calls scoring functions and stores the values in the approriate column.

    Parameters
    ----------
    reading_df : pd.DataFrame
        The reading dataframe
    model_df : pd.DataFrame
        The model dataframe
    graph : nx.DiGraph
        directed graph of the model, necessary for calling kind_score module
    counter: dict
        A dictionary for counting the corrobrated and contradicted interaction
        defulat value is None
    kind_values : dict
        Dictionary assigning Kind Score values
        Default values found in kind_dict
    match_values : dict
        Dictionary assigning Match Score values
        Default values found in match_dict
    attributes : list
        List of attributes compared between the model and the machine reading output
        Default is None
    classify_scheme: str
        The scheme of the classification
        Default value is '1'
    Returns
    -------
    scored = reading_df : pd.DataFrame
        reading dataframe with added scores
    """

    #Create new DF columns for score calculations
    scored_reading_df = reading_df.copy()
    scored_reading_df['Match Score'] = pd.Series()
    scored_reading_df['Kind Score'] = pd.Series()
    scored_reading_df['Epistemic Value'] = pd.Series()
    scored_reading_df['Total Score'] = pd.Series()
    print(reading_df.shape[0])
    #Calculate scores
    for x in range(reading_df.shape[0]):
        scored_reading_df.at[x,'Match Score'] = match_score(x,reading_df,model_df,embedding_match, match_values)
        scored_reading_df.at[x,'Kind Score'] = kind_score(x,model_df,reading_df,graph,embedding_match, counter,kind_values,attributes,classify_scheme,mi_cxn)
        scored_reading_df.at[x,'Epistemic Value'] = epistemic_value(x,reading_df)
        scored_reading_df.at[x,'Total Score'] =  ((scored_reading_df.at[x,'Evidence Score']*scored_reading_df.at[x,'Match Score'])+scored_reading_df.at[x,'Kind Score'])*scored_reading_df.at[x,'Epistemic Value']

    return scored_reading_df
