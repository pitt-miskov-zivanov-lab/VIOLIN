"""
scoring.py

Handles the Match Score, Kind Score, and Epistemic Value functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd
from ast import literal_eval

from numeric import find_element, compare
from network import path_finding

kind_dict = {"strong corroboration" : 2, 
                "weak corroboration1" : 1,
                "weak corroboration2" : 1,
                "weak corroboration3" : 1,
                "hanging extension" : 40, 
                "full extension" : 40, 
                "internal extension" : 40, 
                "specification" : 30, 
                "dir contradiction" : 10,
                "sign contradiction" : 10,
                "att contradiction" : 10,
                "flagged1" : 20,
                "flagged2" : 20,
                "flagged3" : 20}
match_dict = {"source present" : 1, 
                "target present" : 100, 
                "both present" : 10, 
                "neither present" : 0.1}
atts_list = []

def match_score(x,reading_df,model_df,reading_cols,match_values = match_dict):
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
    reading_cols : dict
        Column Header names taken on input
    match_values : dict
        Dictionary assigning Match Score values
        Default values found in match_dict

    Returns
    -------
    match : int
        Match Score value
    """

    #presence of regulated/regulator elements defaults to False
    regulated = False
    regulator = False
    if reading_df.loc[x,reading_cols['pos_source_name']]=="nan": reg_sign = 'neg'
    else: reg_sign = 'pos'

    #Search for regulated from reading in model
    if (find_element("name",reading_df.loc[x,reading_cols['target_name']],reading_df.loc[x,reading_cols['target_type']],model_df) != -1 or 
        find_element("id",reading_df.loc[x,reading_cols['target_id']],reading_df.loc[x,reading_cols['target_type']],model_df) != -1):
        regulated = True

    #Search for regulator from reading in model
    if (find_element("name",reading_df.loc[x,reading_cols[reg_sign+'_source_name']],reading_df.loc[x,reading_cols[reg_sign+'_source_type']],model_df) != -1 or 
        find_element("id",reading_df.loc[x,reading_cols[reg_sign+'_source_id']],reading_df.loc[x,reading_cols[reg_sign+'_source_type']],model_df) != -1):
        regulator = True

    #Scoring Definition
    #only the regulator in the model
    if not regulated and regulator: match = match_values['source present']
    #only the regulated in the model
    elif regulated and not regulator: match = match_values['target present']
    #both regulator and regulated in the model
    elif regulated and regulator: match = match_values['both present']
    #neither in the model
    else: match = match_values['neither present']

    return match


def kind_score(x, model_df, reading_df, graph, reading_cols,
               kind_values = kind_dict, attributes = atts_list, mi_cxn = 'd'):
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
    reading_cols : dict
        Column Header names taken on input
    kind_values : dict
        Dictionary assigning Kind Score values
        Default values found in kind_dict
    attributes : list
        List of attributes compared between the model and the machine reading output
        Default is None
    mi_cxn : str
        What connection type should be assigned to model interactions if not available
        Accepted values are "d" (direct) or "i" (indirect)
        Deafult is "d"

    Returns
    -------
    kind : int
        Kind Score score value
    """

    ###Finding LEE attributes###
    #Finding LEE regulator sign
    if len(list(set(['Reg Sign','Regulator Sign']) & set(reading_df.columns.tolist()))) > 0:
        reg_col = list(set(['Reg Sign','Regulator Sign']) & set(reading_df.columns.tolist()))[0]
        signs = ['Negative','Positive']
        if reading_df.at[x,reg_col][0].lower() in ['activate','positive','increase']: reg_sign = 'Positive'
        else: reg_sign = 'Negative'
        signs.remove(reg_sign)
        opp_sign = signs[0]
    else:
        if reading_df.at[x,reading_cols['pos_source_name']] == "nan": 
            reg_sign = 'Negative'
            opp_sign = 'Positive'
        else:
            reg_sign = 'Positive'
            opp_sign = 'Negative'
    #Finding LEE Connection Type (if not in LEE input, default to indirect, "i")
    if 'Connection Type' in reading_df.columns: lee_cxn_type = reading_df.at[x,'Connection Type']
    else: lee_cxn_type = "i"

    #Finding LEE Other Attributes:
    if len(attributes) > 0:
        #Attribute Column Header Names
        att_cols = [reading_cols['target_'+z.lower().replace(" ","_")] for z in attributes]+[reading_cols[reg_sign[0:3].lower()+'_source_'+z.lower().replace(" ","_")] for z in attributes]
        #Attributes from LEE index 'x'
        reading_atts = [reading_df.at[x,att] for att in att_cols]
    else: reading_atts = []

    ###Comparing to model###
    #Both regulator (source) and regulated (target) node found in the model
    if ((find_element("name",reading_df.at[x,reading_cols['target_name']],reading_df.at[x,reading_cols['target_type']],model_df)!=-1 or 
        find_element("id",reading_df.at[x,reading_cols['target_id']],reading_df.at[x,reading_cols['target_type']],model_df)!=-1) and 
        (find_element("name",reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_name']],reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_type']],model_df)!=-1 or 
        find_element("id",reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_id']],reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_type']],model_df)!=-1)):

        
        #Find indices of regulated element (target) in model 
        if (find_element("name",reading_df.at[x,reading_cols['target_name']],reading_df.at[x,reading_cols['target_type']],model_df)!=-1): 
            model_t_indices = find_element("name",reading_df.at[x,reading_cols['target_name']],reading_df.at[x,reading_cols['target_type']],model_df)
        else: model_t_indices = find_element("id",reading_df.at[x,reading_cols['target_id']],reading_df.at[x,reading_cols['target_type']],model_df)

        #Find indices of regulator element (source) in model
        if (find_element("name",reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_name']],reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_type']],model_df)!=-1): 
            model_s_indices = find_element("name",reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_name']],reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_type']],model_df)
        else: model_s_indices = find_element("id",reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_id']],reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_type']],model_df)
        # print(t_indices,s_indices)

        #Convert regulators into variable names (for path finding)
        model_s_vars = [model_df.loc[i,'Variable'] for i in model_s_indices]

        kinds = []
        #Have to loop over each instance of the regulated element (target) and regulator element (source) in the model, 
        # because the same element may exist in multiple "versions" (usually different locations)
        for t_idx in model_t_indices:
            #Regulator list in model
            model_s_list = model_df.loc[t_idx,reg_sign+" Regulators"]
            #Regulator list of opposite sign
            model_s_opp = model_df.loc[t_idx,opp_sign+" Regulators"]

            for s_idx in model_s_indices:

                # MI with Matched direction, Matched Sign
                if str(model_s_list) != "nan" and model_df.loc[s_idx,'Variable'] in literal_eval(model_s_list):
                    #Index of regulator name within regulator list
                    s_index = literal_eval(model_s_list).index(model_df.loc[s_idx,'Variable'])
                    #Finding index MI regulator variable
                    model_s_variable = literal_eval(model_df.loc[t_idx,reg_sign+' Regulators'])[s_index]
                    model_s_element = list(model_df['Variable']).index(model_s_variable)

                    #Find MI connection type
                    if (reg_sign+' Connection Type') in model_df.columns.values.tolist() and \
                            str(model_df.loc[t_idx,reg_sign+' Connection Type']).lower() not in ['', 'nan', 'NaN']:
                        #Connection type
                        mi_cxn_type = model_df.loc[t_idx,reg_sign+' Connection Type'].split(",")[s_index]
                    else: mi_cxn_type = mi_cxn

                    #List of model attributes to compare to reading attributes
                    model_atts = [model_df.at[t_idx,att] for att in attributes]+[model_df.at[model_s_element,att] for att in attributes]

                    #If LEE ="I" and MI = "I" or LEE = "D" and MI = "D": check attributes
                    if (lee_cxn_type == "i" and mi_cxn_type == "i") or (lee_cxn_type == "d" and mi_cxn_type != "i"):
                        
                        compare_atts = compare(model_atts, reading_atts)
                        #Strong Corroboration - perfect match
                        if compare_atts == 0: kinds.append(kind_values['strong corroboration'])
                        #Specification - the LEE presents new information
                        elif compare_atts == 1: kinds.append(kind_values['specification'])
                        #Weak corroboration - the LEE presents less information than the model interaction
                        elif compare_atts == 2: kinds.append(kind_values['weak corroboration1'])
                        #Contradiction - the LEE presents information that disputes the model interaction
                        elif compare_atts == 3: kinds.append(kind_values['att contradiction'])

                    #If LEE = "D" and MI = "I"
                    elif lee_cxn_type == "d" and mi_cxn_type == "i":
                        compare_atts = compare(model_atts, reading_atts)
                        #If attributes are non-contradictory: LEE is a specification
                        if compare_atts in [0,1,2]: kinds.append(kind_values['specification'])
                        #Else: LEE is a contradiction
                        elif compare_atts == 3: kinds.append(kind_values['att contradiction'])

                    #If LEE ="I" and MI = "D":
                    elif lee_cxn_type == "i" and mi_cxn_type == "d":
                        compare_atts = compare(model_atts, reading_atts)
                        #If attributes are non-contradictory: LEE is a weak corroboration
                        if compare_atts in [0,1,2]: kinds.append(kind_values['weak corroboration2'])
                        #Else: LEE is a contradiction
                        elif compare_atts == 3: kinds.append(kind_values['att contradiction'])

                # MI with Matched direction, Mismatched sign
                elif str(model_s_opp) != "nan" and model_df.loc[s_idx,'Variable'] in literal_eval(model_s_opp) != 0:
                    #LEE is a Sign Contradiction, regardless of connection type
                    kinds.append(kind_values['sign contradiction'])

                # MI with Mismatched direction, Matched sign
                elif (model_df.loc[s_idx,reg_sign+" Regulators"] != "nan" and model_df.loc[t_idx,'Variable'] in literal_eval(model_df.loc[s_idx,reg_sign+" Regulators"])):
                    reg_index = literal_eval(model_df.loc[s_idx,reg_sign+" Regulators"]).index(model_df.loc[t_idx,'Variable'])
                    #Finding connection type
                    if (reg_sign+' Connection Type') in model_df.columns.values.tolist() and \
                            str(model_df.loc[s_idx,reg_sign+' Connection Type']).lower() not in ['', 'nan', 'NaN']:
                        #Connection type
                        mi_cxn_type = model_df.loc[s_idx,reg_sign+' Connection Type'].split(",")[reg_index]
                    else: mi_cxn_type = mi_cxn
                    
                    #Finding index MI regulator variable
                    model_reg_variable = literal_eval(model_df.loc[s_idx,reg_sign+' Regulators'])[reg_index]
                    model_reg_element = list(model_df['Variable']).index(model_reg_variable)
                    #List of model attributes to compare to reading attributes
                    model_atts = [model_df.at[s_idx,att] for att in attributes]+[model_df.at[model_reg_element,att] for att in attributes]

                    # LEE = "I" and MI = "I"
                    if lee_cxn_type == "i" and mi_cxn_type == "i": kinds.append(kind_values['dir contradiction'])
                    
                    # LEE = "D" and MI = "D"
                    elif lee_cxn_type == "d" and mi_cxn_type != "i": 
                        compare_atts = compare(model_atts, reading_atts)
                        #If the attributes are not contradictory - Flagged for manual review
                        if compare_atts in [0,1,2]: kinds.append(kind_values['flagged1'])
                        #Else - Contradiction
                        else: kinds.append(kind_values['dir contradiction'])
                    
                    # LEE = "I" and MI = "D"
                    elif lee_cxn_type == "i" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        #If the attributes are not contradictory - Flagged for manual review
                        if compare_atts in [0,1,2]: kinds.append(kind_values['flagged1'])
                        #Else - Contradiction
                        else: kinds.append(kind_values['dir contradiction'])
                    
                    # LEE = "D" and MI = "I"
                    elif lee_cxn_type == "d" and mi_cxn_type == "i": kinds.append(kind_values['dir contradiction'])
          
                #MI with Mismatched direction, Mismatched sign
                elif (model_df.loc[s_idx,opp_sign+" Regulators"] != "nan" and model_df.loc[t_idx,'Variable'] in literal_eval(model_df.loc[s_idx,opp_sign+" Regulators"])):
                    reg_index = literal_eval(model_df.loc[s_idx,opp_sign+" Regulators"]).index(model_df.loc[t_idx,'Variable'])
                    #Finding connection type
                    if (opp_sign+' Connection Type') in model_df.columns.values.tolist()and \
                            str(model_df.loc[s_idx,reg_sign+' Connection Type']).lower() not in ['', 'nan', 'NaN']:
                        mi_cxn_type = model_df.loc[s_idx,opp_sign+' Connection Type'].split(",")[reg_index]
                    else: mi_cxn_type = mi_cxn

                    #Finding index MI regulator variable
                    model_reg_variable = literal_eval(model_df.loc[s_idx,opp_sign+' Regulators'])[reg_index]
                    model_reg_element = list(model_df['Variable']).index(model_reg_variable)

                    #List of model attributes to compare to reading attributes
                    model_atts = [model_df.at[s_idx,att] for att in attributes]+[model_df.at[model_reg_element,att] for att in attributes]

                    # LEE = "D" and MI = "D"
                    if lee_cxn_type == "d" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        #If the attributes are not contradictory - Flagged for manual review
                        if compare_atts in [0,1,2]: kinds.append(kind_values['flagged1'])
                        #Else - Contradiction
                        else: kinds.append(kind_values['dir contradiction'])
                    # LEE = "D" and MI = "D"
                    elif lee_cxn_type == "d" and mi_cxn_type == "i": kinds.append(kind_values['dir contradiction'])
                    # LEE = "D" and MI = "D"
                    elif lee_cxn_type == "i" and mi_cxn_type != "i":
                        compare_atts = compare(model_atts, reading_atts)
                        #If the attributes are not contradictory - Flagged for manual review
                        if compare_atts in [0,1,2]: kinds.append(kind_values['flagged1'])
                        #Else - Contradiction
                        else: kinds.append(kind_values['dir contradiction'])
                    # LEE = "D" and MI = "D"
                    elif lee_cxn_type == "i" and mi_cxn_type == "i": kinds.append(kind_values['dir contradiction'])
                                    
                #If model does not contain interaction - check for path
                else: 
                    kinds.append(path_finding(model_df.loc[s_idx,'Variable'],model_df.loc[t_idx,'Variable'],reg_sign,model_df,graph,kind_values,lee_cxn_type,reading_atts,attributes))

        if len(kinds) == 1: kind = kinds[0]
        #Strong Corroboration
        elif kind_values['strong corroboration'] in kinds: kind = kind_values['strong corroboration']
        #Weak Corroboration
        elif kind_values['weak corroboration1'] in kinds: kind = kind_values['weak corroboration1']
        elif kind_values['weak corroboration2'] in kinds: kind = kind_values['weak corroboration2']
        elif kind_values['weak corroboration3'] in kinds: kind = kind_values['weak corroboration3']
        #Contradiction
        elif kind_values['dir contradiction'] in kinds: kind = kind_values['dir contradiction']
        elif kind_values['sign contradiction'] in kinds: kind = kind_values['sign contradiction']
        elif kind_values['att contradiction'] in kinds: kind = kind_values['att contradiction']
        #Extensions
        elif kind_values['hanging extension'] in kinds: kind = kind_values['hanging extension']
        elif kind_values['internal extension'] in kinds: kind = kind_values['internal extension']
        elif kind_values['full extension'] in kinds: kind = kind_values['full extension']
        elif kind_values['specification'] in kinds: kind = kind_values['specification']
        #Flagged
        elif kind_values['flagged1'] in kinds: kind = kind_values['flagged1']
        elif kind_values['flagged2'] in kinds: kind = kind_values['flagged2']
        else: kind = kind_values['flagged3']


    #Both Extension - Both nodes from reading not in model
    elif ((find_element("name",reading_df.at[x,reading_cols['target_name']],reading_df.at[x,reading_cols['target_type']],model_df)==-1 and 
        find_element("id",reading_df.at[x,reading_cols['target_id']],reading_df.at[x,reading_cols['target_type']],model_df)==-1) and 
        (find_element("name",reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_name']],reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_type']],model_df)==-1 and 
        find_element("id",reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_id']],reading_df.at[x,reading_cols[reg_sign[0:3].lower()+'_source_type']],model_df)==-1)):
        kind = kind_values['full extension']
    #Hanging Extension - One from reading not in model   
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


def score_reading(reading_df, model_df, graph, reading_cols,
                  kind_values = kind_dict, match_values = match_dict,
                  attributes = atts_list, mi_cxn = 'd'):
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
    reading_cols : dict
        Column Header names taken upon input
    kind_values : dict
        Dictionary assigning Kind Score values
        Default values found in kind_dict
    match_values : dict
        Dictionary assigning Match Score values
        Default values found in match_dict
    attributes : list
        List of attributes compared between the model and the machine reading output
        Default is None

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
        scored_reading_df.at[x,'Match Score'] = match_score(x,reading_df,model_df,reading_cols,match_values)
        scored_reading_df.at[x,'Kind Score'] = kind_score(x,model_df,reading_df,graph,reading_cols,kind_values,attributes,mi_cxn)
        scored_reading_df.at[x,'Epistemic Value'] = epistemic_value(x,reading_df)
        scored_reading_df.at[x,'Total Score'] =  ((scored_reading_df.at[x,'Evidence Score']*scored_reading_df.at[x,'Match Score'])+scored_reading_df.at[x,'Kind Score'])*scored_reading_df.at[x,'Epistemic Value']

    return scored_reading_df
