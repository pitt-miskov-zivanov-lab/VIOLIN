
import nltk
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
    else:
        raise NotImplementedError(f'Metric {metric} is not implemented yet.')
    
