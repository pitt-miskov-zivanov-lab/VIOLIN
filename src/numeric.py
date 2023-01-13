"""
numeric.py

Handles the element finding and comparison functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd
 
def find_element(search_type, element_name, element_type, model_df):
    """
    This function finds the correct indices of an element within the model.
    Because elements can exists as multiple types (Protein, RNA, gene, etc.),
    this function checks the element name/ID along with the element type.
    Function may return a list, if a given element of a specific type 
    exists with varying attributes (such as different locations)

    Parameters
    ----------
    search_type : str
        Whether the element is being searched for by 'name' or 'ID'
    element_name: str
        The name (or ID) of the element being searched for
    element_type: str
        The type of element searched for ('protein', 'protein family', etc.)
    model_df: pd.DataFrame
        The model dataframe

    Returns
    -------
    location : list
        All rows of the model spreadsheet in which the element is found (returns -1 if not found)
    """

    #Searching for element by name
    if search_type == "name":
        #indices of all instances of an element in the model
        indices = [i for i, x in enumerate(list(model_df['Element Name'])) if element_name in x]

    #Searching for element by ID
    else:
        #indices of all instances of an element in the model
        indices = [i for i, x in enumerate(list(model_df['Element IDs'])) if element_name in x]

    #searching for matching element type
    for idx in indices:
        #if type matches, keep that index
        if (model_df.loc[idx,"Element Type"] == element_type or 
            model_df.loc[idx,"Element Type"] in element_type or 
            element_type in model_df.loc[idx,"Element Type"]): next
        #else if type does not match, remove that index
        else: indices.remove(idx)
    
    #If element has been found, return a list of its locations within the model
    if len(indices) > 0: return indices
    #Value -1: means element not found
    else: return -1


def compare(model_atts,reading_atts):
    """
    Compares a list of model attributes to the corresponding LEE attributes, returns numeric value\n
        Attributes are the same (strong corroboration): 0 \n
        Some or all LEE attributes are missing (weak corroboration): 1\n
        Some or all of the model attributes are missing (specification): 2\n
        One or more model attribue differs from the LEE attributes (contradiction): 3

    Parameters
    ----------
    model_att : list
        List of attributes from the model interaction
    reading_att : list
        List of attributes from the literature extracted event (LEE)

    Returns
    -------
    value : int
        Numerical representation of comparison outcome
    """

    #First, indivudally compare each attribute between the machine reading output and model
    #Add outcome to list
    compare_atts = []
    for model in model_atts:
        for reading in reading_atts:
            #Attributes are the same, or neither are available
            if model == reading: compare_atts.append(0)
            elif (model == "nan") and (reading == "nan"): compare_atts.append(0)
            #Reading attribute is not available
            elif (model != "nan") and (reading == "nan"): compare_atts.append(1)
            #Model attribute is not available
            elif (model == "nan") and (reading != "nan"): compare_atts.append(2)
            #Else: both model and reading attributes are available, but they differ
            else: compare_atts.append(3)

    #Final outcome determined by the set of comparisons in compare_atts[]
    #Strong Corroboration - perfect match over all attributes
    if 0 in compare_atts and len(set(compare_atts)) == 1: value = 0

    #Weak Corroboration - model contains more information than the reading for all attributes
    elif 1 in compare_atts and len(set(compare_atts)) == 1: value = 1
    #Weak Corroboration - attributes either match or the model contains more information than the LEE
    elif 0 in compare_atts and 1 in compare_atts and len(set(compare_atts)) == 2: value = 1

    #Specification - LEE contains one or more attributes the model does not
    elif 2 in compare_atts and len(set(compare_atts)) == 1: value = 2
    #Specification - attributes either match or the LEE contains additional attributes
    elif 0 in compare_atts and 2 in compare_atts and len(set(compare_atts)) == 2: value = 2
    #Specification - model or LEE contains more information, depending on attribute
    elif 1 in compare_atts and 2 in compare_atts and len(set(compare_atts)) == 2: value = 2
    
    #Contradiction - some or all of the attributes differ
    elif 3 in compare_atts and len(set(compare_atts)) == 1: value = 3

    #Specification - combination of perfect match, model contains more information, reading contains more information
    else: value = 2

    return value
