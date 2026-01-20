"""
visualize_violin.py

Creates visual representation of VIOLIN output
Updated by April 2025
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
from typing import List

from violin.scoring import KIND_DICT_A, KIND_DICT_B, MATCH_DICT

COLOR_CONFIG = {
    'corroboration': ['#235490','#2B65AD','#718EC5','#B0BEDA', '#E7EFF9'],
    'extension': ['#24552D','#49A155','#7DBA84'],
    'contradiction': ['#BE9735','#E2B441','#F9E6A9'],
    'flagged': ['#AA5626','#C7652E','#F1C6A8', '#F5D4BE', '#FCF1EA'],
}


LABEL_CONFIG = {
    'corroboration': ['strong corroboration', 'empty attribute', 'indirect interaction', 'path corroboration', 'specification'],
    'extension': ['full extension', 'hanging extension', 'internal extension'],
    'contradiction': ['dir contradiction', 'sign contradiction', 'att contradiction'],
    'flagged': ['dir mismatch', 'path mismatch', 'self-regulation', 'flagged4', 'flagged5'],
}


_LABEL_COLORS = {
    'strong corroboration': '#235490',
    'empty attribute': '#2B65AD',
    'indirect interaction': '#718EC5',
    'path corroboration': '#B0BEDA',
    'specification': '#E7EFF9',
    'full extension': '#24552D',
    'hanging extension': '#49A155',
    'internal extension': '#7DBA84',
    'dir contradiction': '#BE9735',
    'sign contradiction': '#E2B441',
    'att contradiction': '#F9E6A9',
    'dir mismatch': '#AA5626',
    'path mismatch': '#C7652E',
    'self-regulation': '#F1C6A8',
    'flagged4': '#F5D4BE',
    'flagged5': '#FCF1EA',
}

_LABEL_DISPLAY = {
    'strong corroboration': "strong corroboration",
    'empty attribute': "empty attribute",
    'indirect interaction': "indirect interaction",
    'path corroboration': "path corroboration",
    'specification': "specification",
    'full extension': "full extension",
    'hanging extension': "hanging extension",
    'internal extension': "internal extension",
    'dir contradiction': "direction contradiction",
    'sign contradiction': "sign contradiction",
    'att contradiction': "attribute contradiction",
    'dir mismatch': "direction mismatch",
    'path mismatch': "path mismatch",
    'self-regulation': "self-regulation",
    'flagged4': "flagged4",
    'flagged5': "flagged5",
}



class ViolinPlot:

    """
    This creates figures of the VIOLIN output:
    evidence score, match score, and total score,
    and classification breakdown

    Parameters
    ----------
    match_values : dict
        Dictionary assigning Match Score Values.
    kind_values : dict
        Dictionary assigning Kind Score values.
    file_name : string
        VIOLIN output to be visualized. Can be specific classification,
        or choosing 'TotalOutput' file will visualize all VIOLIN output.
    filter_opt : str
        How much VIOLIN output should be visualized. Can be filtered
        by top % of total score, evidence score (Se) threshold, or
        total score (St) threshold
        Accepted options are 'X%','Se>Y', or 'St>Z',
        where X, Y, and Z, are values.
        Default is '100%' (Total Output).
    """

    def __init__(self, 
            file_name: str,
            filter_opt: str='100%',
            match_values: dict=None,
            kind_values: dict=None,
            classify_scheme: str='1'):
        
        self.filename = file_name
        self.filter_opt = filter_opt
        assert (classify_scheme in ['1', '2', '3'])
        if match_values is None:
            self.match_values = MATCH_DICT

        if kind_values is None:
            if classify_scheme in ['1', '2']:
                self.kind_values = KIND_DICT_A
            elif classify_scheme == '3':
                self.kind_values = KIND_DICT_B
            else: pass

        # Input file
        self.output = pd.read_csv(file_name, sep=',',index_col=None).fillna("nan")
        if len(self.output) == 0:
            print("The file to plot is empty.")
            return None

        # Filtering by %
        if '%' in filter_opt:
            filter_value = int(filter_opt.replace('%',''))/100
            percent = int(self.output.shape[0]*filter_value)
            self.kept = self.output.head(percent)

        # Filtering by Total Score
        elif 'St>' in filter_opt:
            filter_value = int(filter_opt.replace('St>',''))
            self.kept = self.output.loc[(self.output['Total Score'] >= filter_value)]

        # Filtering by Evidence Score
        elif 'Se>' in filter_opt:
            filter_value = int(filter_opt.replace('Se>',''))
            self.kept = self.output.loc[(self.output['Evidence Score'] >= filter_value)]

        # Else - filter value error
        else:
            raise ValueError('Filter value not accepted'+'\n'+
            'Accepted options are \'X%\',\'Se>Y\', or \'St>Z\','+'\n'+
            'where X, Y, and Z, are numerical values')
        
        self.static_colors = ['royalblue','limegreen','gold','darkorange']
        self.cat_labs = ['corroboration', 'extension', 'contradiction', 'flagged']

    def get_pie_plots(self) -> None:
        """
        This creates figures of the VIOLIN output: the classification distribution shown in pie charts
        """
        # Check if the table is empty
        # cat_vals = []
        # for each in ['dir contradiction','sign contradiction','att contradiction']: cat_vals += [self.kind_values[each]]
        # if len(set(cat_vals))>=1:
        if not self.kept.empty:
            kind = list(self.kept['Kind Score'])
            # Initialize the figure 
            plt.figure(figsize=(16, 4))
            # Create subplots for each category
            for i, cat in enumerate(self.cat_labs):
                plt.subplot(1, 4, i + 1)               
                captions, count_cat = self.count_subcategory(LABEL_CONFIG[cat])
                if all(c == 0 for c in count_cat):
                    plt.pie([1], colors=['grey'])
                    plt.legend(labels=[f"No {cat}"], bbox_to_anchor=((i+1)*0.2, 0), loc="lower center",
                            bbox_transform=plt.gcf().transFigure)
                else:
                    plt.pie(count_cat, colors=COLOR_CONFIG[cat])
                    plt.legend(labels=captions, bbox_to_anchor=((i+1)*0.2, 0), loc="lower center",
                            bbox_transform=plt.gcf().transFigure)
        else:
            raise ValueError("The file to plot is empty.")

    def get_summary_plots(self, save=False, merge=True) -> None:
        """
        A summary plot composed with category distribution, evidence score, match score, and total score
        """

        # Split the kept dataframe into categories
        cat_dfs = self.get_category_df()

        # Get overall plot of the data 
        plt.figure(figsize=(12, 6))
        X_axis = np.arange(len(LABEL_CONFIG))
        cat_count = {}
        cat_count = {cat: cat_df.shape[0] for cat, cat_df in cat_dfs.items()}
        # bar plot for the main category distribution
        plt.subplot(2, 2, 1)
        plt.bar(X_axis, list(cat_count.values()), label=list(cat_count.keys()), color=self.static_colors)
        plt.xticks(X_axis, list(cat_count.keys()))
        plt.ylabel('Number of IOLs')

        # Evidence Score plots
        scores = list(set(self.kept['Evidence Score']))
        X_axis = np.arange(len(scores))
        
        assert ("Evidence Score" in self.kept.columns)
        plt.subplot(2, 2, 2)
        self.get_category_score_plot(list(self.kind_values.keys()), "Evidence Score", "Evidence Score", merge_subcat=merge, save=False)

        # Match Score plots
        assert ("Match Score" in self.kept.columns)
        plt.subplot(2, 2, 3)
        self.get_category_score_plot(list(self.kind_values.keys()), "Match Score", "Match Score", merge_subcat=merge, save=False)

        # Total Score plots
        assert ("Total Score" in self.kept.columns)
        plt.subplot(2, 2, 4)
        self.get_category_score_plot(list(self.kind_values.keys()), "Total Score", "Total Score", merge_subcat=merge, save=False)

        plt.tight_layout()
        if save:
            plt.savefig('Output_Overview.png',bbox_inches = "tight",dpi=200)
        plt.close()

    def get_category_summary(self, category: str, save_name: str='', save=False) -> None:
        """
        Plot the score (evidence, match, total) for specified categories.
        """
        if category not in LABEL_CONFIG:
            raise ValueError(f"Category {category} not found in {LABEL_CONFIG.keys()}.")

        assert ("Evidence Score" in self.kept.columns)
        sub_cat_list = set(LABEL_CONFIG[category]).intersection(self.kind_values.keys())
        plt.subplot(2, 2, 1)
        self.get_category_score_plot(list(sub_cat_list), "Evidence Score", "Evidence Score", save=False)

        # Match Score plots
        assert ("Match Score" in self.kept.columns)
        plt.subplot(2, 2, 2)
        self.get_category_score_plot(list(sub_cat_list), "Match Score", "Match Score", save=False)

        # Total Score plots
        assert ("Total Score" in self.kept.columns)
        plt.subplot(2, 2, 3)
        self.get_category_score_plot(LABEL_CONFIG[category], "Total Score", "Total Score", save=False)

        # category pie plot
        plt.subplot(2, 2, 4)
        captions, count_cat = self.count_subcategory(LABEL_CONFIG[category])
        if all(c == 0 for c in count_cat):
            plt.pie([1], colors=['grey'])
            plt.legend(labels=[f"No {category}"], bbox_to_anchor=(1, 0), loc="lower right",
                    bbox_transform=plt.gcf().transFigure)
        else:
            plt.pie(count_cat, colors=COLOR_CONFIG[category])
            plt.legend(labels=captions, bbox_to_anchor=(1, 0), loc="lower right",
                    bbox_transform=plt.gcf().transFigure)
        
        plt.tight_layout()
        if save:
            save_name = category if save_name == '' else save_name
            plt.savefig(f'{save_name}.png',bbox_inches = "tight",dpi=200)
        plt.close()
        return 
    
    # --------------------------------------Internal utility functions-------------------------------------- #

    def get_category_score_plot(self, sub_cat_list: List[str], score_name:str, save_name:str='', save=False, merge_subcat:bool=False) -> None:
        """
        Plot the score (evidence, match, total) for specified categories.
        """
        if any(sub_cat not in _LABEL_DISPLAY.keys() for sub_cat in sub_cat_list):
            raise ValueError(f"Category {sub_cat_list} not found in {_LABEL_DISPLAY.keys()}.")

        # Get the category dataframe
        dfs = self.get_category_df(sub_cats=sub_cat_list)  

        # Merge subcategories if specified
        if merge_subcat:
            cat_dfs = {}
            for cat in self.cat_labs:
                dfs_to_concat = [
                    dfs[sub_cat]
                    for sub_cat in sub_cat_list
                    if sub_cat in LABEL_CONFIG[cat] and not dfs[sub_cat].empty
                ]
                if dfs_to_concat:  # Only concat if list is not empty
                    cat_dfs[cat] = pd.concat(dfs_to_concat, ignore_index=True)
                else:
                    cat_dfs[cat] = pd.DataFrame()
        else:
            cat_dfs = dfs

        scores = list(set(self.kept[score_name]))
        X_axis = np.arange(len(scores)) 
        for idx, (name, cat_df) in enumerate(cat_dfs.items()):

            if merge_subcat:
                color = self.static_colors[self.cat_labs.index(name)]
                label = name
            else:
                color = _LABEL_COLORS[name]
                label = _LABEL_DISPLAY[name]

            # Get the subcategory dataframe
            counts = self._count_score_for_category(self.kept, cat_df, score_name)
            plt.bar(X_axis+(idx*0.2),counts,0.2,color=color,label=label)
        step = max(1, len(scores) // 10)
        plt.xticks(X_axis[::step], sorted([scores[i] for i in range(0, len(scores), step)]), rotation=45)
        plt.yscale('log')
        plt.legend(prop={'size': 6}) 
        plt.ylabel('Number of IISs')
        plt.xlabel(score_name)
        plt.tight_layout()
        if save:
            plt.savefig(f'{save_name}.png',bbox_inches = "tight",dpi=200)
        
    def count_subcategory(self, sub_cat_list: List[str]) -> int:
        if any(x not in _LABEL_DISPLAY.keys() for x in sub_cat_list):
            raise ValueError(f"Category {sub_cat_list} not a subset of {_LABEL_DISPLAY.keys()}.")
        
        kind = self.kept['Kind Score'].tolist()
        sub_cat_count = {}; sub_cat_plot_caption = []
        for sub_cat in sub_cat_list:
            # Skip flagged 4 and 5 if not in kind_values
            if sub_cat not in list(self.kind_values.keys()):
                pass
            else:
                # Get count for each sub-category
                sub_cat_count[sub_cat] = kind.count(self.kind_values[sub_cat])
                # Get legend and label for each sub-category
                sub_cat_plot_caption.append(f"{_LABEL_DISPLAY[sub_cat]}: {sub_cat_count[sub_cat]}")
        
        return sub_cat_plot_caption, np.array(list(sub_cat_count.values()))
    
    def get_category_df(self, sub_cats:List[str]=None) -> dict:
        cat_df_dict = {}
        if sub_cats is None:
            for cat, sub_cats in LABEL_CONFIG.items():
                cat_df_dict[cat] = self.kept[self.kept['Kind Score'].isin([self.kind_values[sub_cat] for sub_cat in sub_cats\
                                                                           if sub_cat in self.kind_values.keys()])]
                
        else:
            for sub_cat in sub_cats:
                if sub_cat in self.kind_values.keys():
                    cat_df_dict[sub_cat] = self.kept[self.kept['Kind Score'].isin([self.kind_values[sub_cat]])]
        return cat_df_dict
    
    @staticmethod
    def _count_score_for_category(df:pd.DataFrame, category_df: pd.DataFrame, score_name: str) -> List[int]:
        """
        Count scores (evidence, match, total) for each category
        """
        score = category_df[score_name].value_counts().keys().tolist()
        counts = category_df[score_name].value_counts().tolist()
        for value in list(df[score_name]):
            # Add placeholder if value not in score
            if value not in score:
                score += [value]
                counts += [0]
        counts = [x for _,x in sorted(zip(score, counts))]
        return counts