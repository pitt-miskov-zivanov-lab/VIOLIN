# This is a script for summarizing the results from all different readers
# Usage:
# python summarization.py --result_dir ./output --out_dir ./ --reader_name GPT --scheme v1

import pandas as pd
import glob
import re
import sys
import os
import argparse

sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), os.pardir, 'src/violin/')))
from violin.in_out import preprocessing_model, KIND_DICT_A, KIND_DICT_B

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

# FILES = [
#     "translated_SkeMel133_biorecipe",
#     "translated_ModelB_discrete_biorecipe",
#     "translated_SkeMel133_biorecipe_combined_10contradictions",
#     "translated_ModelB_discrete_biorecipe_combined_10contradictions",
#     "translated_SkeMel133_biorecipe_combined_10randoms",
#     "translated_ModelB_discrete_biorecipe_combined_10randoms",
# ]


def merge_duplicates(reading_df, col_names):
    """
    This function merges duplicate interactions and calculates evidence score of each LEE

    Parameters
    ----------
    reading_df : pd.DataFrame
        The dataframe of the machine reading output
    col_names: list
        Specifically the column headings used to determine if interactions are identical

    Returns
    -------
    counted_reading : pd.DataFrame
        A new dataframe with the evidence count and PMCID list for each interaction
    """

    #Convert reading to lower case, to prevent issues with case difference
    reading = reading_df.apply(lambda x: x.astype(str).str.lower())
    #The columns that aren't used to determine duplicates (such as Paper ID or Evidence Text)
    remainder = [x for x in reading_df.columns if x not in col_names]

    #As VIOLIN Identifies duplicates, it merges attributes from the remainder list into a single cell
    #This is how we count the number of times an IOL appears, and keep track of paper IDs and evidence text
    counted_reading = reading.groupby(col_names)[remainder[0]].apply(list).reset_index(name=remainder[0])
    for x in range(1,len(remainder)):
        sub = reading.groupby(col_names)[remainder[x]].apply(list).reset_index(name=remainder[x])
        counted_reading[remainder[x]] = sub[remainder[x]]

    return counted_reading


def main(result_dir: str, out_dir: str, reader_name: str, scheme: str ) -> None:
        
        global FILES
        
        print("Make sure your results match the scheme you are used for summary...")
        print(f'current result folder: {result_dir}')
        print(f'scheme: {scheme}')

        count = {}

        if 'GPT' in reader_name or 'LLAMA' in reader_name:
            FILES = FILES[1:]

        for f in FILES:
            print(f)

            for c_ in ['corroboration', 'contradiction', 'flagged', 'extension']:
                file = os.path.join(result_dir, f'{f}_{c_}.csv')
                out_name = f

                df = pd.read_csv(file, index_col=None).fillna('nan').astype(str)
                df = df.map(lambda x: x.lower().strip() if isinstance(x, str) else x)

                if c_ not in count:
                    count[c_] = {}

                merged_df = merge_duplicates(df,
                                             ["Regulator Name",
                                              "Regulator Type",
                                              "Regulator ID",
                                              "Regulated Name",
                                              "Regulated Type",
                                              "Regulated ID",
                                              "Sign"
                ])

                direct_cnx = len([x for row in range(len(merged_df))
                                        for x in merged_df.loc[row, 'Connection Type'] if x.lower() == 'd'])
                mechanism_attr = len([x for row in range(len(merged_df))
                                        for x in merged_df.loc[row, 'Mechanism'] if x.lower() not in ['none', 'nan']])
                phos_attr = len([x for row in range(len(merged_df))
                                        for x in merged_df.loc[row, 'Mechanism'] if x.lower() == 'phosphorylation'])

                if 'connection' not in count[c_]: count[c_]['connection'] = {}
                if 'mechanism' not in count[c_]: count[c_]['mechanism'] = {}
                if 'phosphorylation' not in count[c_]: count[c_]['phosphorylation'] = {}

                count[c_]['connection'][out_name] = direct_cnx
                count[c_]['mechanism'][out_name] = mechanism_attr
                count[c_]['phosphorylation'][out_name] = phos_attr

                if scheme in ['v1', 'v2']:
                    kind_dict = KIND_DICT_A
                else:
                    kind_dict = KIND_DICT_B

                if c_ == 'corroboration':
                    if 'strong corroboration' not in count[c_]:
                        count[c_]['strong corroboration'] = {}
                    if 'empty attribute' not in count[c_]:
                        count[c_]['empty attribute'] = {}
                    if 'indirect interaction' not in count[c_]:
                        count[c_]['indirect interaction'] = {}
                    if 'path corroboration' not in count[c_]:
                        count[c_]['path corroboration'] = {}
                    if 'specification' not in count[c_]:
                        count[c_]['specification'] = {}

                    count[c_]['strong corroboration'][out_name] = len(df[df['Kind Score'] == str(kind_dict['strong corroboration'])])
                    count[c_]['empty attribute'][out_name] = len(df[df['Kind Score'] == str(kind_dict['empty attribute'])])
                    count[c_]['indirect interaction'][out_name] = len(df[df['Kind Score'] == str(kind_dict['indirect interaction'])])
                    count[c_]['path corroboration'][out_name] = len(df[df['Kind Score'] == str(kind_dict['path corroboration'])])
                    count[c_]['specification'][out_name] = len(df[df['Kind Score'] == str(kind_dict['specification'])])

                elif c_ == 'contradiction':
                    if "sign contradiction" not in count[c_]:
                        count[c_]["sign contradiction"] = {}
                    if 'dir contradiction' not in count[c_]:
                        count[c_]['dir contradiction'] = {}
                    if 'att contradiction' not in count[c_]:
                        count[c_]['att contradiction'] = {}
                    # Placeholders for redundent columns
                    if 'contradiction4' not in count[c_]:
                        count[c_]['contradiction4'] = {}
                    if 'contradiction5' not in count[c_]:
                        count[c_]['contradiction5'] = {}

                    count[c_]["sign contradiction"][out_name] = len(df[df['Kind Score'] == str(kind_dict["sign contradiction"])])
                    count[c_]['dir contradiction'][out_name] = len(df[df['Kind Score'] == str(kind_dict['dir contradiction'])])
                    count[c_]['att contradiction'][out_name] = len(df[df['Kind Score'] == str(kind_dict['att contradiction'])])
                    count[c_]['contradiction4'][out_name] = 0
                    count[c_]['contradiction5'][out_name] = 0

                elif c_ == 'flagged':
                    if "dir mismatch" not in count[c_]:
                        count[c_]["dir mismatch"] = {}
                    if "path mismatch" not in count[c_]:
                        count[c_]["path mismatch"] = {}
                    if "self-regulation" not in count[c_]:
                        count[c_]["self-regulation"] = {}
                    if "flagged4" not in count[c_]:
                        count[c_]["flagged4"] = {}
                    if 'flagged5' not in count[c_]:
                        count[c_]["flagged5"] = {}

                    count[c_]["dir mismatch"][out_name] = len(df[df['Kind Score'] == str(kind_dict["dir mismatch"])])
                    count[c_]['path mismatch'][out_name] = len(df[df['Kind Score'] == str(kind_dict['path mismatch'])])
                    count[c_]['self-regulation'][out_name] = len(df[df['Kind Score'] == str(kind_dict['self-regulation'])])
                    if scheme in ['v1', 'v2']:
                        count[c_]['flagged4'][out_name] = 0
                        count[c_]['flagged5'][out_name] = 0
                    else:
                        count[c_]['flagged4'][out_name] = len(df[df['Kind Score'] == str(kind_dict['flagged4'])])
                        count[c_]['flagged5'][out_name] = len(df[df['Kind Score'] == str(kind_dict['flagged5'])])

                else:
                    if "full extension" not in count[c_]:
                        count[c_]["full extension"] = {}
                    if "hanging extension" not in count[c_]:
                        count[c_]["hanging extension"] = {}
                    if "internal extension" not in count[c_]:
                        count[c_]["internal extension"] = {}
                    # Placeholders for redundent columns
                    if "extension4" not in count[c_]:
                        count[c_]["extension4"] = {}
                    if 'extension5' not in count[c_]:
                        count[c_]["extension5"] = {}

                    count[c_]["full extension"][out_name] = len(df[df['Kind Score'] == str(kind_dict["full extension"])])
                    count[c_]['hanging extension'][out_name] = len(df[df['Kind Score'] == str(kind_dict['hanging extension'])])
                    count[c_]['internal extension'][out_name] = len(df[df['Kind Score'] == str(kind_dict['internal extension'])])
                    count[c_]['extension4'][out_name] = 0
                    count[c_]['extension5'][out_name] = 0


        for _ in ['corroboration', 'contradiction', 'flagged', 'extension']:

            dict_ = {}
            #dict_['reading'] = TEST_FILES
            dict_['reading'] = FILES
            for key, value in count[_].items():
                if key == 'reading':
                    pass
                else:
                    dict_[key] = list(count[_][key].values())

            count_df = pd.DataFrame(dict_)
            out_file = os.path.join(out_dir, reader_name, f'{scheme}_{_}_summary.csv')
            os.makedirs(os.path.dirname(out_file), exist_ok=True)
            count_df.to_csv(out_file, index=False)


if __name__ == '__main__':
    args = argparse.ArgumentParser(description="Classification results.")
    args.add_argument('--result_dir', required=True, help="The classified results folder name.")
    args.add_argument('--reader_name', required=True, help="The name of the reader.")
    args.add_argument('--out_dir', required=True, help="The directory of output summary.")
    args.add_argument('--scheme', required=True, help="The classification scheme, i.e. v1, v2, or v3")
    args = args.parse_args()
    
    main(args.result_dir, args.out_dir, args.reader_name, args.scheme)
