"""
scoring.py

Handles the Match Score, Kind Score, and Epistemic Value functions for VIOLIN
Created November 2019 - Casey Hansen MeLoDy Lab
"""

import pandas as pd
import networkx as nx
from violin.formatting import get_type
from violin.numeric import find_element
from violin.decision_tree import kind_score
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
                include_hits=False,
                name_match_threshold: float=None,
                match_sim_metric: str='jaccard',
                normalize_type: bool=False,
                differ_gene_type: bool=False) -> pd.DataFrame:
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
    name_match_threshold : float
        The confidence threshold for finding element in the model.
    match_sim_metric : str
        The similarity metric for calculating match score.
    normalize_type: bool
        Whether to normalize the type of elements in the model and reading dataframes before scoring.
    differ_gene_type: bool
        Whether to distinguish gene type and protein type 

    Returns
    -------
    scored = reading_df : pd.DataFrame
        reading dataframe with added scores.
    """
    assert (classify_scheme in ['1', '2', '3'])
    assert (mi_cxn in ['d', 'i'])

    logger.info(f"metric: {match_sim_metric}")
    logger.info(f"name_match_threshold: {name_match_threshold}")

    if kind_values is None:
        if classify_scheme in ['1', '2']:
            kind_values = KIND_DICT_A

        elif classify_scheme == '3':
            kind_values = KIND_DICT_B
    else:
        pass

    if match_values is None:
        match_values = MATCH_DICT

    # To compatible with UI, we also add a normalization step during classification
    # Normalize the type of elements in the model and reading dataframes
    if normalize_type: 
        reading_df['Regulator Type'] = reading_df['Regulator Type'].apply(lambda x: get_type(x, differ_gene=differ_gene_type))
        reading_df['Regulated Type'] = reading_df['Regulated Type'].apply(lambda x: get_type(x, differ_gene=differ_gene_type))
        model_df['Element Type'] = model_df['Element Type'].apply(lambda x: get_type(x, differ_gene=differ_gene_type))


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
        scored_reading_df.at[x,'Kind Score'], interaction_hit = kind_score(
            x,model_df,
            reading_df,
            graph, 
            counter,
            kind_values,
            attributes,
            classify_scheme,
            mi_cxn, 
            name_match_threshold=name_match_threshold,
            match_sim_metric=match_sim_metric)
        scored_reading_df.at[x,'Epistemic Value'] = epistemic_value(x,reading_df)
        scored_reading_df.at[x,'Total Score'] =  ((scored_reading_df.at[x,'Evidence Score']*scored_reading_df.at[x,'Match Score'])+int(scored_reading_df.at[x,'Kind Score']))*scored_reading_df.at[x,'Epistemic Value']

        if interaction_hit != None:
            collections.append({"category": interaction_hit[0], 
                            "kind_score": str(interaction_hit[1]),
                            "target_index": interaction_hit[2],
                            "source_index": interaction_hit[3],
                            "interaction_index": x,
                            "regulated_conf": interaction_hit[4],
                            "regulator_conf": interaction_hit[5]})
            prob = interaction_hit[4] * interaction_hit[5]
        
        else:
            prob = -1 # Unknown error

        scored_reading_df.at[x, 'Entity Match Score'] = round(prob, 2)

    col = scored_reading_df.pop("Entity Match Score")
    scored_reading_df.insert(0, "Entity Match Score", col)
    if include_hits:
        return scored_reading_df, collections
    else: 
        return scored_reading_df
