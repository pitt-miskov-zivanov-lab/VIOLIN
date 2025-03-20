import os
import sys
import re
import argparse
import glob
import json 
from tqdm import tqdm
import time

import requests
import httplib2 as http
try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

import pandas as pd
import numpy as np
import gilda

from violin.formatting import get_subtype_abbr
from collections import Counter

headers = {
    "Accept": 'application/json'
}

h = http.Http()

# model columns that are able to be normalized
model_columns = ['Element Name', 'Cell Line', 'Cell Type', 'Tissue Type', 'Organism']

BioRECIPE_reading_col = ["Regulator Name", "Regulator Type", "Regulator Subtype", "Regulator HGNC Symbol",
                         "Regulator Database", "Regulator ID", "Regulator Compartment", "Regulator Compartment ID",
                         "Regulated Name", "Regulated Type", "Regulated Subtype", "Regulated HGNC Symbol",
                         "Regulated Database", "Regulated ID", "Regulated Compartment", "Regulated Compartment ID",
                         "Sign", "Connection Type", "Mechanism", "Site",
                         "Cell Line", "Cell Type", "Tissue Type", "Organism",
                         "Score", "Source", "Statements", "Paper IDs"]

type_abbr_dict = {
             'proteinfamily': 'pf',
             'proteincomplex': 'pf',
             'protein': 'pn',
             'chemical': 'che',
             'chemicalfamily': 'cf',
             'gene': 'gene',
             'rna': 'rna',
             'mutation': 'mut',
             'biologicalprocess': 'bp',
             'complex': 'complex',
             'nan': 'nan'
}




class ModelProofreading:
    @staticmethod
    def string2list(x: str) -> list:
        if x.lower() in ['', 'nan', 'none']:
            return []
        
        else:
            x = re.split(r'[/,]\s*', x)
            return x
        
    def __init__(self, model_filename: str):
        self.model = pd.read_excel(model_filename).fillna("nan")
        self.listname = self._get_listname()
        self.repeat_ele = None
        self.pos_not_found = []
        self.neg_not_found = []
        self.row_not_unique = []
        self.row_cxn_corr = []

    def _get_listname(self) -> list:
        # Define the "Listname" for every element in model file
        # Forward synthesize the name in regulator list
        listname = []
        name = list(self.model['Element Name'])
        etype = [type_abbr_dict[x].strip() for x in self.model['Element Type'].fillna('nan'). str.lower().str.replace(' ', '')]
        subtype = [get_subtype_abbr(x) for x in self.model['Element Subtype'].fillna('nan')]
        location = [x.replace(':', '').strip() for x in self.model['Compartment ID'].fillna('nan')]
        listname = [f'{a}_{b}_{c}_{d}' for a, b ,c ,d, e in zip(name, etype, subtype, location, range(len(name)))]
        return listname
        
        
    def proofreading(self) -> None:
        """
        A function to check if all the regulators in the columns are correct:
            - check the existence of the regulators 
            - check the uniqueness of the elements in the model file
            - check the uniqueness of regulators for every element 
            - check the correspondence of connection type of regulators
        """
        # Check if the elements are unique or not
        listname_freq = Counter(self.listname)
        
        assert all( x >= 1 for x in listname_freq.values())

        if len(listname_freq) != self.listname:
            self.repeat_ele = [listname for listname, freq in list(listname_freq.items()) if freq != 1]
            raise ValueError(f"(1/4)the elements in the file are not unique, existing repeating nodes: {",".join(self.repeat_ele)}\n")
        else:
            pass

        # Check the existence of the regulators

        for row in range(len(self.model)):
            for x in self.string2list(self.model.loc[row, "Positive Regulator List"]):
                x = x.strip()
                if x in self.listname:
                    pass 
                else: 
                    self.pos_not_found.append('row {row} -> {x}(pos)')
            
            for x in self.string2list(self.model.loc[row, "Negative Regulator List"]):
                x = x.strip()
                if x in self.listname:
                    pass
                else:
                    self.neg_not_found.append('row {row} -> {y}(neg)')
        
        if len(self.pos_not_found) or len(self.neg_not_found): 
            
            print(f"(2/4)regulators not found: {",".join(self.pos_not_found)}, {",".join(self.neg_not_found)}\n")
        else:
            pass

        # Check the uniqueness of the regulators 
        for row in range(len(self.model)):
            pos_list = self.string2list(self.model.loc[row, "Positive Regulator List"])
            neg_list = self.string2list(self.model.loc[row, "Negative Regulator List"])
            if len(pos_list) != len(set(pos_list)) or len(pos_list) != len(set(neg_list)):
                self.row_not_unique.append(row)
        
        if self.row_not_unique:
            print(f"(3/4)Regulators are not unique in row: {",".join(self.row_not_unique)}\n")
        else:
            pass

        # Check the correspondence of the reuglators
        for row in range(len(self.model)):
            pos_list = self.string2list(self.model.loc[row, "Positive Regulator List"])
            neg_list = self.string2list(self.model.loc[row, "Negative Regulator List"])
            pos_cxn_list = self.string2list(self.model.loc[row, "Positive Connection Type List"])
            neg_cxn_list = self.string2list(self.model.loc[row, "Negative Connection Type List"])
            if len(pos_list) != len(pos_cxn_list) or len(neg_list) != len(neg_cxn_list):
                print(f"(4/4)Regulator are not unique in row: {",".join(self.row_not_unique)}\n")


def get_hgnc_symbol(hgnc_id: str, hgnc_dict: dict={}, url: str='https://rest.genenames.org/fetch/hgnc_id') -> str:
    """
    A function to convert HGNC id to HGNC symbol
    
    Parameters
    ----------
    hgnc_id: str
        An HGNC id string that adheres to the pattern specified at https://bioregistry.io/registry/hgnc.
    hgnc_dict: dict
        A cached dict store the HGNC symbol and id
    url: str
        the url for HGNC's rest api.
    
    Returns
    -------
    symbol: str
        A HGNC symbol corresponds to the HGNC id.
    """

    # Check if it is a local unique identifier
    if re.match(r'^\d{1,5}$', hgnc_id):
        pass
    
    # Check if it is a CURIE 
    elif re.match(r'^hgnc:\d{1,5}$', hgnc_id):
        hgnc_id = hgnc_id.replace('hgnc:', '')
    
    # Return if is is unsearchable
    else:
        print(f'unable to proceed {hgnc_id}')
        return hgnc_id

    # Check if we can find the number
    if hgnc_id in hgnc_dict.keys():
        symbol = hgnc_dict[hgnc_id]
    else:
        response, content = h.request(
            url + f'/{hgnc_id}',
            'GET',
            '',
            headers
        )
        data = json.loads(content)
        status_code = False; i = 0; symbol = ''
        while status_code != True and i < 10:
            try:
                i += 1
                response, content = h.request(
                    url + f'/{hgnc_id}',
                    'GET',
                    '',
                    headers)
    
                if response['status'] == '200':
                    symbol = data['response']['docs'][0]['symbol']
                    hgnc_dict[hgnc_id] = symbol
                    status_code = True
                else:
                    pass
            except Exception as e:
                print(e)
                time.sleep(1)

    return symbol


def change_id(x:str) -> str:
    """
    Exceptional database identifier switch

    Parameters
    ----------
    x: str
        the identifer

    Returns 
    -------
    x: str
        A processed identifier
    """
    if str(x) == 'nan':
        return x
    else:
        if type(x) is str:
            if ',' in x:
                id_list = []
                for id_x in x.split(','):
                    id_x = change_id(id_x)
                    id_list.append(id_x)
                return ','.join(id_list)
            else:
                if re.match(r'\s*\{*\s*chebi:.*', x.lower()) is not None:
                    return x.lower().replace('chebi:', '')

                elif re.match(r'\s*\{*\s*go:.*', x.lower()) is not None:
                    return x.lower().replace('go:', '')

                elif re.match(r'\s*\{*\s*hmdb.*', x.lower()) is not None:
                    return x.lower().replace('hmdb', '')
                else:
                    return x
        else:
            return x


def normalize_element(df: pd.DataFrame, clns: list=['Regulated', 'Regulator']) -> pd.DataFrame:
    """
    A function of using GILDA,
    Parameters
    ----------
    df: pd.DataFrame
        A dict of input file
    clns: list
        A list of columns of interest

    Returns
    -------
    df: pandas dataframe
        A dataframe of input file where the columns of interest are normalized as CURIEs
    """
    public_url = 'http://grounding.indra.bio.ground'
    with open('hgnc_symbol_dict.json', 'w') as f:
        cached_hgnc_dict = json.load(f)
        
    # Register clns
    for cln in tqdm(clns):
        #gilda.ground_df(df, source_column=cln, target_column=f"Normalized {cln}")
        for row in range(len(df)):
            try:
                res = requests.post(public_url, json={'text':df.loc[row, cln], 'context': df.loc[row, "Statements"]})
                result_list = res.json()
                scores = np.array([x['score'] for x in result_list])
                if result_list:
                    result = result_list[np.argmax(scores)]
                    df.loc[row, f'Normalized {cln} Name'] = result['term']['entry_name']
                    df.loc[row, f'Normalized {cln} ID'] = change_id(result['term']['id'])
                    df.loc[row, f'Normalized {cln} Database'] = result['term']['db']

                    hgnc_id = df.loc[row, '{cln} HGNC Symbol']
                    df.loc[row, f'Normalized {cln} HGNC Symbol'] = get_hgnc_symbol(hgnc_id=hgnc_id, hgnc_dict=cached_hgnc_dict)
                else:
                    pass
            except Exception as e:
                print(e)
                pass 


    return df


def mapping_name_reading(df: pd.DataFrame) -> pd.DataFrame:
    """
    This is a function to normalize the bio-entities in BioRECIPE interaction list file,
    GILDA has a threshold for normalization, we only overwrite the term when confidence is over the threshold.
    Parameters
    ----------
    df: pandas dataframe
        a dataframe for the reading file.
    Returns
    -------
    df: pandas dataframe
        A dict of normalized reading file
    """
    clns = df.columns
    for cln in clns:
        if 'Normalized' in cln:
            org_cln = cln.split('Normalized')[-1].strip()
            df[cln] = df[cln].where(df[cln].isna(), df[cln].str.replace(':', ''))
            df[org_cln] = df[cln].where(df[cln].notna(), df[org_cln])
            df = df.drop(columns=[cln])
        else:
            pass
    return df



def main():
    parser = argparse.ArgumentParser(description="A script for name normalization")

    parser.add_argument("--format", type=str, choices=['model', 'reading'], help="Choose a file type from options.")
    parser.add_argument("--file", type=str, help="The pathname of excel file in BioRECIPE format.", default=None)
    parser.add_argument("--folder", type=str, help="A directory of excel files in BioRECIPE format.", default=None)
    parser.add_argument("--selected_columns", type=str, choices=['Regulator', 'Regulated'], help="please use commas to separate the column names of "
                                                             "interest", default='Regulated, Regulator')
    parser.add_argument("--outdir", type=str, help="A directory to store the outputs.", default=None)

    # Parse the arguments
    args = parser.parse_args()
    # get column names
    columns = [x.strip() for x in args.selected_columns.split(',')]


    # Check if output directory is correct
    if args.outdir is None or (not os.path.isdir(args.outdir)):
        raise ValueError('error in out directory.')

    if args.file is None and args.folder is None:
        raise ValueError("assign a value in file or folder.")

    # Check pathname of the file
    elif args.file is not None:
        # Check if output directory is same as input directory
        if os.path.dirname(args.file) == args.outdir:
            choice = input(
                "program might overwrite the original file, would you like to continue? (y/n)").strip().lower()
            if choice == 'y':
                pass
            else:
                sys.exit()
        else:
            pass

        normalized_df = normalize_element(df=pd.read_excel(args.file, index_col=None), clns=columns)
        f_name = args.file.split('/')[-1]

        if args.format == 'reading':
            normalized_df = mapping_name_reading(normalized_df)

        elif args.format == 'model':
            print(f+'-'*20)
            checker = ModelProofreading(normalized_df)
            checker.proofreading()

        # Save out file
        normalized_df.to_excel(f'{args.outdir}/{f_name}', index=False)
        print(f'finish normalizing {args.file}', end='\n')

    elif args.folder is not None:
        if os.path.isdir(args.folder):
            # Check if output directory is same as input directory
            if args.folder == args.outdir:
                choice = input(
                    "program might overwrite the original file, would you like to continue? (y/n)").strip().lower()
                if choice == 'y':
                    pass
                else:
                    sys.exit()
            else:
                pass

            fs = glob.glob(f'{args.folder}/*.xlsx')
            for f in fs:
                f_df = pd.read_excel(f, index_col=None)
                f_name = f.split('/')[-1]
                normalized_df = normalize_element(df=f_df, clns=columns)
                if args.format == 'reading':
                    normalized_df = mapping_name_reading(normalized_df)
                    normalized_df.to_excel(f'{args.outdir}/{f_name}', index=False)
                elif args.format == 'model':
                    print(f+'-'*20)
                    checker = ModelProofreading(normalized_df)
                    checker.proofreading()

                print(f'finish normalizing {f}', end='\n')
        else:
            raise ValueError('folder is not a directory')

    else:
        raise ValueError('choose either file or folder as input.')


if __name__ == "__main__":
    main()
