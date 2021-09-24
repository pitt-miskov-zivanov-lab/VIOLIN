"""
network.py

Handles the functions necessary for finding paths within the model for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd
import numpy as np
import networkx as nx
from VIOLIN.numeric import compare


def node_edge_list(model_df):
    """
    This function converts the model from the BioRECIPES format into a node-edge list for use with NetworkX

    Parameters
    ----------
    model_df : pd.DataFrame
        The model dataframe, must be in BioRECIPES format

    Returns
    -------
    node_edge_list : nx.DiGraph
        A directed graph representation of the model
    """

    #If elements are defined by variables, use the variable names. Else, use the common names
    if 'Variable' in model_df.columns: target = 'Variable'
    else: target = 'Element Name'

    #Subset of the model, just element and regulator columns
    graph = model_df[[target,'Positive Regulators','Negative Regulators']]
    #removes 'nan' placeholder
    graph = graph.replace('nan','')

    #remove excess punctuation from the regulator cells
    graph['Positive Regulators'] = graph['Positive Regulators'].str.replace('[','').str.replace(']','').str.replace('\'','')
    graph['Negative Regulators'] = graph['Negative Regulators'].str.replace('[','').str.replace(']','').str.replace('\'','')
    
    #combine regulators into one column, separated by '-' symbol
    #positiveRegulators_negativeRegulators
    graph['Regulators'] = graph['Positive Regulators']+'-'+graph['Negative Regulators']
    graph = graph.drop(columns=['Positive Regulators','Negative Regulators'])

    #Split each row by '-' character
    #This allows us to assign weights so that we know the "sign" (positive/negative) of the regulator
    #Even rows are positive regulators, odd rows are negative regulators
    graph = graph.set_index([target]).stack().str.split('-', expand=True).stack().unstack(-2).reset_index(-1, drop=True).reset_index()
    #Assign 'weights' to edge defined edge, 0 for positive regulators, 1 for negative regulators
    for y in range(graph.shape[0]):
        if y%2 == 0: graph.at[y,'weight'] = 0
        else: graph.at[y,'weight'] = 1
    
    #Remove rows without a regulator node (housekeeping)
    graph = graph.replace(r'^\s*$', np.nan, regex=True).dropna()
    graph = graph.set_index([target,'weight']).stack().str.split(',', expand=True).stack().unstack(-2).reset_index(-1, drop=True).reset_index()
    #Remove any lingering whitespace
    graph['Regulators'] = graph['Regulators'].str.strip()
    graph[target] = graph[target].str.strip()
    #Output edgelist
    # graph.to_csv(r'/NodeEdgeList.csv')
    #Create NetworkX directed graph
    node_edge_list = nx.from_pandas_edgelist(graph,'Regulators',target,'weight',create_using=nx.DiGraph())
    return node_edge_list


def path_finding(regulator,regulated,sign,model_df,graph,kind_values,reading_cxn_type,reading_atts,attributes):
    """
    This function searches for a path between the reading regulator and regulated in the model,
    and calculates the kind score based on the results

    Parameters
    ----------
    regulator : str
        Element variable name of the regulator node
    regulated : str
        Element variable name of the regulated node
    sign : str
        Sign of regulated node
    model_df : pd.DataFrame
        Model dataframe
    graph : nx.DiGraph
        Model edgelist to create network for finding paths between elements
    kind_values : dict
        Dictionary containing the numerical values for the Kind Score classifications
    reading_cxn_type : str
        Connection Type of interaction from reading - 'i' for indirect, 'd' for direct

    Returns
    -------
    kind : int
        Kind Score value for the interaction
    """

    #Sign of regulator; assigned same numbering as node_edge_list() function
    #negativeRegulators = 1; positiveRegulators = 0
    if sign == 'Negative': sign = 1
    else: sign = 0
    

    #Have to make sure regulator and regulated are in the directed graph representation of the model
    #Some nodes may be in the model, but aren't regulated/regulators anywhere
    if regulator in graph and regulated in graph:
        #If there is a path of the same direction and LEE = D: internal extension
        if nx.has_path(graph,regulator,regulated) and len(nx.shortest_path(graph,source=regulator,target=regulated))>1 and reading_cxn_type == "d": kind = kind_values['internal extension']
        #If there is a path of the same direction and LEE = I: check sign and attributes
        elif nx.has_path(graph,regulator,regulated) and len(nx.shortest_path(graph,source=regulator,target=regulated))>1 and reading_cxn_type == "i": 
            #Finding atts of beginning and end of path
            s_idx = list(model_df['Variable']).index(regulator)
            t_idx = list(model_df['Variable']).index(regulated)
            model_atts = [model_df.at[t_idx,att] for att in attributes]+[model_df.at[s_idx,att] for att in attributes]
            compare_atts = compare(model_atts, reading_atts)

            #Finding Path sign
            #path list
            path = nx.shortest_path(graph,source=regulator,target=regulated,weight='weight')
            #Check path sign
            path_wgt = 0
            idx = 0
            #Sum the edge weights to determine the overall effect
            while idx < len(path)-1:
                path_wgt += graph[path[idx]][path[idx+1]]['weight']
                idx += 1
            # if %2 = 0, then positive regulation, if %2 = 1, then negative regulation
            # Weak corroboration - regulation matches reading 
            if path_wgt%2 == sign and compare_atts in [0,1,2]: kind = kind_values['weak corroboration3']
            #Flagged - Regulation same sign, but contradictory attributes
            elif path_wgt%2 == sign and compare_atts == 3: kind = kind_values['flagged2']
            # Flagged - Regulation opposite sign as reading
            else: kind = kind_values['flagged2']

        #If there is a path of the opposite direction - Flagged
        elif nx.has_path(graph,regulated,regulator) and len(nx.shortest_path(graph,source=regulated,target=regulator))>1:
            kind = kind_values['flagged2']
        
        #If there is a self regulation (regulator is both target and source)
        elif nx.has_path(graph,regulator,regulated) and len(nx.shortest_path(graph,source=regulator,target=regulated))==1: 
            kind = kind_values['flagged3']

        #If there is no path
        else: kind = kind_values['internal extension']
    else: kind = kind_values['internal extension']

    return kind
