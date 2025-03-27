import pandas as pd
import glob
import os
import argparse
import time

from violin.numeric import find_element
from violin.in_out import preprocessing_model, preprocessing_reading, output, KIND_DICT_A, KIND_DICT_B
from violin.network import node_edge_list
from violin.scoring import score_reading, MATCH_DICT


FILES_TEST = ["RB1"]

FILES = ["RA1",
         "RA2",
         "RA3",
         "RA4",
         "RB1",
         "RB2",
         "RB3",
         "RB_star_1",
         "RB_star_2"]


evidence_scoring_cols = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol", "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                        "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol", "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                        "Sign", "Connection Type", "Mechanism", "Site",
                        "Cell Line", "Cell Type", "Tissue Type", "Organism"]
# attributes = ['Regulated Compartment ID', 'Regulator Compartment ID']


def main(scheme: str, reader: str, out_dir: str, attributes: list) -> None:

    approach = scheme.split('v')[-1]

    if approach in ['1', '2']:
        kind_dict = KIND_DICT_A
    elif approach == '3':
        kind_dict = KIND_DICT_B

    if attributes == None:
        attributes = []
    else: 
        pass


    model_files = ['input/models/SkMel133_biorecipe.xlsx', 'input/models/ModelB_discrete_biorecipe.xlsx']
    reading_A_files = glob.glob(f'input/interactions/{reader}/RA*.xlsx')
    reading_B_files = glob.glob(f'input/interactions/{reader}/RB*.xlsx')

    model_A_df = preprocessing_model(model_files[0])
    model_B_df = preprocessing_model(model_files[1])

    graph_A = node_edge_list(model_A_df)
    graph_B = node_edge_list(model_B_df)

    reader = f'{reader}/{scheme}'

    if 'FLUTE' in reader:
        reading_A_files = [x for x in reading_A_files if '_BioRECIPE_filtered' in x]
        reading_B_files = [x for x in reading_B_files if '_BioRECIPE_filtered' in x]
    else:
        pass


    if os.path.isdir(os.path.join(out_dir, reader)):
        pass
    else:
        os.makedirs(os.path.join(out_dir, reader))
    
    print(attributes)

    for reading_file in reading_A_files:
        if 'FLUTE' in reader:
            output_file = os.path.join(out_dir, reader, reading_file.split('/')[-1].split('_BioRECIPE_filtered')[0])
        else:
            output_file = os.path.join(out_dir, reader, reading_file.split('/')[-1].split('_reading_BioRECIPE')[0])
        print(output_file)
        time1 = time.time()
        reading_df = preprocessing_reading(reading=reading_file, 
                                        evidence_score_cols=evidence_scoring_cols, 
                                        atts=attributes)
        counter_A = {'corroboration': [], 'contradiction': []}
        time2 = time.time()
        scored = score_reading(reading_df, 
                        model_A_df, 
                        graph_A, 
                        counter=counter_A,
                        kind_values=kind_dict, 
                        match_values=MATCH_DICT, 
                        attributes=attributes, 
                        classify_scheme=approach,
                        )
        print(f'classification time: {time.time() - time2}')
        output(scored, output_file, kind_values=kind_dict)
        print(f'total time: {time.time() - time1}')
        print('corroboration in model: {}'.format(len(set(counter_A['corroboration']))))
        print('contradiction in model: {}'.format(len(set(counter_A['contradiction']))))

    for reading_file in reading_B_files:
        if 'FLUTE' in reader:
            output_file = os.path.join(out_dir, reader, reading_file.split('/')[-1].split('_BioRECIPE_filtered')[0])
        else:
            output_file = os.path.join(out_dir, reader, reading_file.split('/')[-1].split('_reading_BioRECIPE')[0])
        print(output_file)
        time1 = time.time()
        reading_df = preprocessing_reading(reading=reading_file, 
                                        evidence_score_cols=evidence_scoring_cols, 
                                        atts=attributes)
        counter_B = {'corroboration': [], 'contradiction': []}
        time2 = time.time()
        scored = score_reading(reading_df, 
                        model_B_df, 
                        graph_B, 
                        counter=counter_B,
                        kind_values=kind_dict, 
                        match_values=MATCH_DICT, 
                        attributes=attributes, 
                        classify_scheme=approach)
        print(f'classification time: {time.time() - time2}')
        output(scored, output_file, kind_values=kind_dict)
        print(f'total time: {time.time() - time1}')
        print('corroboration in model: {}'.format(len(set(counter_B['corroboration']))))
        print('contradiction in model: {}'.format(len(set(counter_B['contradiction']))))

if __name__ == '__main__':
    args = argparse.ArgumentParser(description='Evaluation Scheme')
    args.add_argument('-s', '--scheme', required=False, choices=['v1', 'v2', 'v3'], default='v1', help='Scheme name, i.e. v1, v2, or v3')
    args.add_argument('-r', '--reader', required=True, choices=['REACH', 'INDRA', 'GPT', 'LLAMA', 'FLUTE/REACH',
                    'FLUTE/INDRA', 'FLUTE/GPT', 'FLUTE/LLAMA'], help='Readers name, e.g. REACH.')
    args.add_argument('-o', '--output', required=True, help='Output directory')
    args.add_argument('-a', '--attributes', required=False, choices=[
        'Regulated Compartment ID', 'Regulator Compartment ID', 'Mechanism', 'Cell Line', 'Cell Type', 'Tissue Type', 'Organism'
    ], 
    default=['Regulated Compartment ID', 'Regulator Compartment ID'],nargs='+', help='Element, influence, and context attributes used to compare.')
    args = args.parse_args()
    print(args.attributes)
    main(args.scheme, args.reader, args.output, args.attributes)


