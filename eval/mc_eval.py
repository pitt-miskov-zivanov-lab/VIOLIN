import os
import sys
import json

import pandas as pd
import numpy as np
import argparse 

from violin.in_out import preprocessing_model, preprocessing_reading, output, BioRECIPE_READING_COL, KIND_DICT_A
from violin.network import node_edge_list
from violin.scoring import score_reading
from violin.visualize_violin import LABEL_CONFIG

EVIDENCE_COLS = ["Score", "Source", "Statements", "Paper IDs"]
EVIDENCE_SCORING_COLS = [x for x in BioRECIPE_READING_COL if x not in EVIDENCE_COLS]

ATTRIBUTES = ["Regulated Compartment ID", "Regulator Compartment ID"]

def parse_int(value):
    if isinstance(value, str):
        return int(value)
    elif isinstance(value, int):
        return value
    elif isinstance(value, float):
        return int(value)
    elif isinstance(value, list):
        return parse_int(value[0])
    else:
        raise ValueError(f"Cannot convert {value} to int")
    
def overlap(a, b):
    """
    Check if two lists have any common elements.
    """
    return any(item in a for item in b)
    

def count_stat(results, scenario='Class'):
    assert scenario in ['Class', 'Kind Score'], f'Unknown scenario {scenario}'
    stat = {}
    stat['attributes'] = {}
    stat['categories'] = {}
    _score2label = {v: k for k, v in KIND_DICT_A.items()}
    results[scenario] = results[scenario].apply(lambda x: parse_int(x) if isinstance(x, (list, str)) else x)
    for cls in LABEL_CONFIG.keys():
        stat['attributes'][cls] = {}
        stat['categories'][cls] = {}
        sub_cls_scores = [KIND_DICT_A[x] for x in LABEL_CONFIG[cls] if x in KIND_DICT_A]
        
        cls_rows = results[results[scenario].isin(sub_cls_scores)]
        stat['attributes'][cls]['Direct connection'] = cls_rows[cls_rows['Connection Type'] == 'd'].shape[0]
        stat['attributes'][cls]['Mechanism'] = cls_rows[cls_rows['Mechanism'] != 'none'].shape[0]
        stat['attributes'][cls]['Phosphorylation'] = cls_rows[cls_rows['Mechanism'] == 'phosphorylation'].shape[0]

        stat['categories'][cls] = {_score2label[subcat]: cls_rows[cls_rows[scenario] == subcat].shape[0] for subcat in sub_cls_scores}
    return stat

def evaluate(results, subcat=False):
    """
    Evaluate the results against the ground truth.
    """
    gt_labels = results['Class']
    predict_labels = results['Kind Score']                                                                               

    metric = {}
    sub_cls_scores = {}

    for cls in LABEL_CONFIG.keys():
        sub_cls = [x for x in LABEL_CONFIG[cls] if x in KIND_DICT_A]
        sub_cls_scores[cls] = [int(KIND_DICT_A[x]) for x in sub_cls]

    for cls, sub_cls_score in sub_cls_scores.items():
        tp, tn, fp, fn = [], [], [], []
        for idx, (gt_label, predict_label) in enumerate(zip(gt_labels, predict_labels)):
            if not isinstance(gt_label, list):
                gt_label = [gt_label]

            gt_label = [parse_int(x) for x in gt_label]
            for scores in sub_cls_scores.values():
                if gt_label[0] in scores:
                    gt_label = scores
                    break

            # True Positives
            if overlap(gt_label, sub_cls_score) and predict_label in sub_cls_score:
                tp.append(1)
            # True Negatives
            elif not overlap(gt_label, sub_cls_score) and predict_label not in sub_cls_score:
                tn.append(1)
            # False Positives
            elif not overlap(gt_label, sub_cls_score) and predict_label in sub_cls_score:
                fp.append(1)
            # False Negatives
            elif overlap(gt_label, sub_cls_score) and predict_label not in sub_cls_score:
                fn.append(1)
            else:
                raise ValueError(f'Unknown case for {cls} at index {idx} with gt_labels {gt_label} and predict_label {predict_label}')

        recall = sum(tp) / (sum(tp) + sum(fn)) if (sum(tp) + sum(fn)) > 0 else 0
        precision = sum(tp) / (sum(tp) + sum(fp)) if (sum(tp) + sum(fp)) > 0 else 0
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

        metric[cls] = {
            'TP': sum(tp),
            'TN': sum(tn),
            'FP': sum(fp),
            'FN': sum(fn),
            'Precision': precision,
            'Recall': recall,
            'F1 Score': f1_score
        }
    
    if subcat:
        metric['Subcategories'] = {}
        gt_labels = gt_labels.apply(lambda x: parse_int(x) if isinstance(x, (list, str)) else x)
        for subcls, score in KIND_DICT_A.items():
            tp=[]; tn=[]; fp=[]; fn=[]
            for row in results.iterrows():
                row = row[1]
                if parse_int(row['Class']) == score and row['Kind Score'] == score:
                    tp.append(1)
                elif parse_int(row['Class']) != score and row['Kind Score'] == score:
                    fp.append(1)
                elif parse_int(row['Class']) == score and row['Kind Score'] != score:
                    fn.append(1)
                else:
                    tn.append(1)

            tp = sum(tp); fp = sum(fp); fn = sum(fn); tn = sum(tn)
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

            metric['Subcategories'][subcls] = {
                'Precision': precision,
                'Recall': recall,
                'F1 Score': f1_score
            }



    return metric

def classify_reading(reading_file, model_file, approach='1', evidence_scoring_cols=EVIDENCE_SCORING_COLS, attributes=ATTRIBUTES):

    model_df = preprocessing_model(model_file)
    reading_df = preprocessing_reading(reading=reading_file, evidence_score_cols=evidence_scoring_cols, atts=attributes)
    graph = node_edge_list(model_df)

    print(f'# of interactions: {len(pd.read_excel(reading_file))}')
    scored = score_reading(reading_df, model_df, graph, attributes=attributes, classify_scheme=approach)
    print(len(scored))
    return scored

def main():
    args = argparse.ArgumentParser(description='Evaluation Scheme')
    args.add_argument('--model', required=True, help='Path to the model file')
    args.add_argument('--reading', required=True, help='Path to the reading file')
    args.add_argument('--output', required=True, help='Path to the output file')

    args = args.parse_args()

    model_file = args.model
    output_file = args.output

    # Predict the results
    results = classify_reading(reading_file=args.reading, model_file=model_file)

    # Evaluate the results
    evaluation_results = evaluate(results, subcat=True)
    # evaluation_results = evaluate_metric(results)
    # Count statistics
    ground_stats = count_stat(results, scenario='Class')
    prediction_stats = count_stat(results, scenario='Kind Score')
    # Summarize the results
    summary = {
        'ground_truth': ground_stats,
        'prediction': prediction_stats,
        'evaluation': evaluation_results
    }
    # Save the classified results
    basename = os.path.basename(args.reading).split('.')[0]
    output(results, f'./examples/extension_manually_checking/curation_correction_250529/output/{basename}', classify_scheme='1')

    # Save the results
    print(evaluation_results)
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=4)

    

if __name__ == '__main__':
    main()


