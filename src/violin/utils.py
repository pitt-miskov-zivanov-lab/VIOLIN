
import nltk
import gilda
from typing import List, Union, Tuple, Dict, Any
from bioregistry import parse_curie
import logging 

logging.basicConfig(level=logging.INFO)
_RESERVED_NULL = ['null', 'none', 'nan', 'n/a', 'na']

def get_similarity_score(name1:str, name2:str, metric: str='jaccard')-> float:
    """
    Get the confidence score of two names,
    currently using edit distance as the metric
    """
    if metric == 'edit_sim':
        total_len = max(len(name1), len(name2))
        # special case: when both names are null-like, we consider them as a meaningless match
        if total_len == 0: 
            return 0.0
        if name1.lower() in _RESERVED_NULL and name2.lower() in _RESERVED_NULL:
            return 0.0
        distance = nltk.edit_distance(name1, name2, transpositions=False)
        score = 1 - distance / total_len
        return score
    elif metric == 'jaccard':
        total_len = max(len(name1), len(name2))
        # special case: when both names are null-like, we consider them as a meaningless match
        if total_len == 0: 
            return 0.0
        if name1.lower() in _RESERVED_NULL and name2.lower() in _RESERVED_NULL:
            return 0.0
        set1 = set(name1); set2 = set(name2)
        return len(set1 & set2) / len(set1 | set2)
    
    elif metric == 'none':
        return 0.0
    else:
        raise NotImplementedError(f'Metric {metric} is not implemented yet.')
    

def gilda_grounding(name: str, 
                    context:str=None,
                    namespaces:List[str]=None,
                    normalize_curie: bool=True) -> Union[Tuple[Dict[str, str], float], Tuple[List[str], float], Tuple[None, None]]: 
    
    """
    wrapper for gilda grounding API, which returns the best match and its score for the name

    """
    data = gilda.ground(
        name, context=context, namespaces=namespaces
    )

    if data is not None and len(data) > 0:
        id_dict = {}
        data = data[0]
        db_ids = data.get_grounding_dict()
        conf = data.score 
        if normalize_curie:
            for db_id in db_ids:
                db, id_ = parse_curie(db_id)
                id_dict[db] = id_
            return id_dict, conf
        else:
            return list(db_ids.keys()), conf
    
    else:
        return None, None
    

def gilda_find_element(iname, icontext=None, inamespace=None, mname_list=None, model_df=None, threshold=0.5) -> Tuple[Any, Any]: 
    # mname_list: a dict including meta information of mname_list {curie: {name: List[str], model_index: List[int], confidence: List[float]}
    threshold = threshold if threshold is not None else 0.99
    if model_df is not None: 
        raise NotImplementedError("Model-based disambiguation is not implemented yet.")
    
    else:
        id_dict, conf = gilda_grounding(
            iname, context=icontext, namespaces=inamespace, normalize_curie=False
        )
        if id_dict is None: 
            return (-1, [])
        
        else: 
            indices = []; confs = []
            for id_ in id_dict:
                if id_ in mname_list: 
                    ele_dict = mname_list[id_]
                    for name, m_idx, m_conf in zip(ele_dict['name'], ele_dict['model_index'], ele_dict['confidence']):
                        c = m_conf * conf
                        if c >= threshold:
                            indices.append(m_idx)
                            confs.append(c)
            if len(indices) == 0: 
                return (-1, [])

            return (indices, confs)
        

def get_model_element_curies(model_df) -> dict:

    model_name_dict = {}
    for i, name in enumerate(model_df['Element Name']):
        ids, conf = gilda_grounding(name, normalize_curie=False)
        if ids is not None:
            for id_ in ids:
                if id_ in model_name_dict:
                    model_name_dict[id_]['name'].append(str(name))
                    model_name_dict[id_]['model_index'].append(i+2)
                    model_name_dict[id_]['confidence'].append(conf)
                else:
                    model_name_dict[id_] = {'name': [str(name)], 'model_index': [i], 'confidence': [conf]}
        elif ids is None: 
            ids = model_df.loc[i, 'Element IDs']
            model_name_dict[ids] = {'name': [str(name)], 'model_index': [i], 'confidence': [1]}
    return model_name_dict

# A class utilize entity grounding to find the matched element
class ElementMatcher: 
    
    def __init__(self, model_df=None): 
        if model_df is not None:
            self.model_df = model_df
            self.model_name_dict = get_model_element_curies(model_df)
        else:
            self.model_df = None
            self.model_name_dict = {}

    def initialize_model_curies(self, model_df):

        if self.model_df is not None:
            logging.info("Model data already initialized, replacing original model data in ElementMatcher.")

        self.model_df = model_df
        self.model_name_dict = get_model_element_curies(model_df)

    def find_element(self, iname, icontext=None, inamespace=None, threshold=0.5) -> Tuple[Any, Any]:
        if self.model_df is None or self.model_name_dict is None:
            raise ValueError("Model data is not initialized.")
        
        return gilda_find_element(iname, 
                                  icontext, 
                                  inamespace, 
                                  mname_list=self.model_name_dict, 
                                  model_df=None, 
                                  threshold=threshold)