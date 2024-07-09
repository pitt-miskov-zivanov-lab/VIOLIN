import pandas as pd
import glob
import re
import sys
import os
from ast import literal_eval

sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), os.pardir, 'src/violin/')))
from violin.numeric import find_element
from violin.in_out import preprocessing_model

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

# FILES = ["RA1",
#          "RA2",
#          "RA3",
#          "RA4",
#          "RB1",
#          "RB2",
#          "RB3",
#          "RB_star_1",
#          "RB_star_2",
#          "RA2_0_1",
#          "RA2_0_1_1",
#          "RB2_1",
#          "RB2_0_1"]

# TEST_FILES = ["translated_ModelB_discrete",
#               "translated_SkeMel133",
#               "translated_SkeMel133_biorecipe_combined_10contradictions",
#               "translated_ModelB_discrete_biorecipe_combined_10contradictions",
#               "translated_SkeMel133_biorecipe_combined_10randoms",
#               "translated_ModelB_discrete_biorecipe_combined_10randoms"]

TEST_FILES = ["RA2",
              "RA21",
              "RA2_0_1",
              "RA2_0_1_1",
              "RB2",
              "RB21",
              "RB2_0_1"
              ]

kind_dict_a = {"strong corroboration" : 2,
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

kind_dict_b = {"strong corroboration" : 2,
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
                "flagged4":17,
                "flagged5":16}


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
    #This is how we count the number of times an LEE appears, and keep track of paper IDs and evidence text
    counted_reading = reading.groupby(col_names)[remainder[0]].apply(list).reset_index(name=remainder[0])
    for x in range(1,len(remainder)):
        sub = reading.groupby(col_names)[remainder[x]].apply(list).reset_index(name=remainder[x])
        counted_reading[remainder[x]] = sub[remainder[x]]

    return counted_reading

if __name__ == '__main__':
    scheme_file = 'v2'
    count = {}
    folder = 'output'
    reader = 'INDRA/v2'
    model_A = preprocessing_model('input/models/SkMel133_biorecipe.xlsx')
    model_B = preprocessing_model('input/models/ModelB_discrete_biorecipe.xlsx')

    filenames = [file for file in glob.glob(scheme_file + '/*.csv') if re.match(r'^((?!score).)*$', file)]
    if 'GPT' in reader or 'LLAMA' in reader:
        FILES = FILES[1:]
    #for f in TEST_FILES:
    for f in FILES:
        print(f)
        if re.search(r'^RA', f) or re.search(r'SkeMel133', f):
            model_df = model_A
        else:
            model_df = model_B

        for c_ in ['corroborations', 'contradictions', 'flagged', 'extensions']:
            file = f'{folder}/{reader}/{f}_{c_}.csv'
            out_name = f

            df = pd.read_csv(file, index_col=None).fillna('nan').astype(str)
            df = df.applymap(lambda x: x.lower().strip() if isinstance(x, str) else x)
            #df['Regulator Name'] = df['Regulator Name'].str.replace('-', '').str.replace(' ', '').str.replace('_', '')

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

            # df = df.drop_duplicates(subset=[
            #     'Regulator Name',
            #     'Regulator Type', 'Regulator Compartment',
            #     'Regulated Name', 'Regulated Compartment',
            #     'Regulated Type',
            #     'Sign', 'Mechanism']
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

            if scheme_file in ['v1', 'v2']:
                kind_dict = kind_dict_a
            else:
                kind_dict = kind_dict_b

            if c_ == 'corroborations':
                if 'strong corr' not in count[c_]:
                    count[c_]['strong corr'] = {}
                if 'weak corr1' not in count[c_]:
                    count[c_]['weak corr1'] = {}
                if 'weak corr2' not in count[c_]:
                    count[c_]['weak corr2'] = {}
                if 'weak corr3' not in count[c_]:
                    count[c_]['weak corr3'] = {}
                if 'Specification' not in count[c_]:
                    count[c_]['Specification'] = {}

                count[c_]['strong corr'][out_name] = len(df[df['Kind Score'] == str(kind_dict['strong corroboration'])])
                count[c_]['weak corr1'][out_name] = len(df[df['Kind Score'] == str(kind_dict['empty attribute'])])
                count[c_]['weak corr2'][out_name] = len(df[df['Kind Score'] == str(kind_dict['indirect interaction'])])
                count[c_]['weak corr3'][out_name] = len(df[df['Kind Score'] == str(kind_dict['path corroboration'])])
                count[c_]['Specification'][out_name] = len(df[df['Kind Score'] == str(kind_dict['specification'])])

            elif c_ == 'contradictions':
                if "sign contradiction" not in count[c_]:
                    count[c_]["sign contradiction"] = {}
                if 'dir contradiction' not in count[c_]:
                    count[c_]['dir contradiction'] = {}
                if 'att contradiction' not in count[c_]:
                    count[c_]['att contradiction'] = {}
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
                count[c_]['flagged4'][out_name] = 0
                count[c_]['flagged5'][out_name] = 0

            else:
                if "full extension" not in count[c_]:
                    count[c_]["full extension"] = {}
                if "hanging extension" not in count[c_]:
                    count[c_]["hanging extension"] = {}
                if "internal extension" not in count[c_]:
                    count[c_]["internal extension"] = {}
                if "extension4" not in count[c_]:
                    count[c_]["extension4"] = {}
                if 'extension5' not in count[c_]:
                    count[c_]["extension5"] = {}

                count[c_]["full extension"][out_name] = len(df[df['Kind Score'] == str(kind_dict["full extension"])])
                count[c_]['hanging extension'][out_name] = len(df[df['Kind Score'] == str(kind_dict['hanging extension'])])
                count[c_]['internal extension'][out_name] = len(df[df['Kind Score'] == str(kind_dict['internal extension'])])
                count[c_]['extension4'][out_name] = 0
                count[c_]['extension5'][out_name] = 0


    for _ in ['corroborations', 'contradictions', 'flagged', 'extensions']:

        dict_ = {}
        #dict_['reading'] = TEST_FILES
        dict_['reading'] = FILES
        for key, value in count[_].items():
            if key == 'reading':
                pass
            else:
                dict_[key] = list(count[_][key].values())

        count_df = pd.DataFrame(dict_)

        count_df.to_csv(f'{folder}/{reader}/{scheme_file}_{_}_summary_id.csv', index=False)


