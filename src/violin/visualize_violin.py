"""
visualize_violin.py

Creates visual representation of VIOLIN output
Created November 2019 - Casey Hansen MeLoDy Lab
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
plt.rcParams["font.family"] = "Arial"

def visualize (match_values, kind_values, file_name, filter_opt='100%'):
    """
    This creates graphs of the VIOLIN output:
    evidence score, match score, and total score,
    and classification breakdown

    Parameters
    ----------
    match_values : dict
        Dictionary assigning Match Score Values
    kind_values : dict
        Dictionary assigning Kind Score values
    file_name : string
        VIOLIN output to be visualized. Can be specific classification,
        or choosing 'TotalOutput' file will visualize all VIOLIN output
    filter_opt : str
        How much VIOLIN output should be visualized. Can be filtered
        by top % of total score, evidence score (Se) threshold, or
        total score (St) threshold
        Accepted options are 'X%','Se>Y', or 'St>Z',
        where X, Y, and Z, are values
        Default is '100%' (Total Output)
    """

    # Input file
    output = pd.read_csv(file_name, sep=',',index_col=None).fillna("nan")

    # Filtering by %
    if '%' in filter_opt:
        filter_value = int(filter_opt.replace('%',''))/100
        percent = int(output.shape[0]*filter_value)
        kept = output.head(percent)

    # Filtering by Total Score
    elif 'St>' in filter_opt:
        filter_value = int(filter_opt.replace('St>',''))
        kept = output.loc[(output['Total Score'] >= filter_value)]

    # Filtering by Evidence Score
    elif 'Se>' in filter_opt:
        filter_value = int(filter_opt.replace('Se>',''))
        kept = output.loc[(output['Evidence Score'] >= filter_value)]

    # Else - filter value error
    else:
        raise ValueError('Filter value not accepted'+'\n'+
        'Accepted options are \'X%\',\'Se>Y\', or \'St>Z\','+'\n'+
        'where X, Y, and Z, are numerical values')

    # If visualizing all categories of output
    if '_outputDF' in file_name:
        # Separating output my category
        corroborations = kept[kept["Kind Score"].isin([kind_values['strong corroboration'],
                                                           kind_values['weak corroboration1'],
                                                           kind_values['weak corroboration2'],
                                                           kind_values['weak corroboration3']])]
        extensions = kept[kept["Kind Score"].isin([kind_values['full extension'],
                                                    kind_values['hanging extension'],
                                                    kind_values['internal extension'],
                                                    kind_values['specification']])]
        contradictions = kept[kept["Kind Score"].isin([kind_values['dir contradiction'],
                                                        kind_values['sign contradiction'],
                                                        kind_values['att contradiction']])]
        flagged = kept[kept["Kind Score"].isin([kind_values['flagged1'],
                                                    kind_values['flagged2'],
                                                    kind_values['flagged3']])]

        plt.figure(figsize=(12, 6))
        # Plot main category distribution
        category = ['Corroboration','Extension','Contradiction','Flagged']
        X_axis = np.arange(len(category))
        mycolors = ['royalblue','limegreen','gold','darkorange']
        counts = [corroborations.shape[0],extensions.shape[0],contradictions.shape[0],flagged.shape[0]]

        plt.subplot(2, 2, 1)
        plt.bar(X_axis,counts,label=category,color=mycolors)
        plt.xticks(X_axis, category)
        plt.ylabel('Number of LEEs')

        #Evidence Score plots
        scores = list(set(kept['Evidence Score']))
        X_axis = np.arange(len(scores))
        category = [corroborations['Evidence Score'],extensions['Evidence Score'],
                    contradictions['Evidence Score'],flagged['Evidence Score']]
        cat_lab = ['corroboration', 'extensions', 'contradictions', 'flagged']
        colors = ['royalblue','limegreen','gold','darkorange']

        plt.subplot(2, 2, 2)
        for idx in range(len(category)):
            evidence = category[idx].value_counts().keys().tolist()
            counts = category[idx].value_counts().tolist()
            # housekeeping to make sure all data is the same size
            for value in scores:
                if value not in evidence:
                    evidence += [value]
                    counts += [0]
            counts = [x for _,x in sorted(zip(evidence,counts))]
            plt.bar(X_axis+(idx*0.2),counts,0.2,color=colors[idx],label=cat_lab[idx])
        plt.xticks(X_axis, scores, rotation=45)
        plt.yscale('log')
        plt.legend(prop={'size': 6})
        plt.ylabel('Number of LEEs')
        plt.xlabel('Evidence Score')

        #Match Score plots
        scores = list(set(kept['Match Score']))
        X_axis = np.arange(len(scores))
        category = [corroborations['Match Score'],extensions['Match Score'],
                    contradictions['Match Score'],flagged['Match Score']]
        cat_lab = ['corroboration', 'extensions', 'contradictions', 'flagged']
        colors = ['royalblue','limegreen','gold','darkorange']

        plt.subplot(2, 2, 3)
        for idx in range(len(category)):
            evidence = category[idx].value_counts().keys().tolist()
            counts = category[idx].value_counts().tolist()
            for value in scores:
                if value not in evidence:
                    evidence += [value]
                    counts += [0]
            counts = [x for _,x in sorted(zip(evidence,counts))]
            plt.bar(X_axis+(idx*0.2),counts,0.2,color=colors[idx],label=cat_lab[idx])
        plt.xticks(X_axis, scores)
        plt.legend(prop={'size': 6})
        plt.ylabel('Number of LEEs')
        plt.xlabel('Match Score')

        #Total Score plots
        scores = list(set(kept['Total Score'].tolist()))
        scores = [float(i.replace('{','').replace('}','')) for i in scores]
        X_axis = np.arange(len(scores))
        category = [corroborations['Total Score'],extensions['Total Score'],
                    contradictions['Total Score'],flagged['Total Score']]
        cat_lab = ['corroboration', 'extensions', 'contradictions', 'flagged']
        colors = ['royalblue','limegreen','gold','darkorange']

        plt.subplot(2, 2, 4)
        for idx in range(4):
            evidence = [float(i.replace('{','').replace('}','')) for i in category[idx].value_counts().keys().tolist()]
            counts = category[idx].value_counts().tolist()
            for value in scores:
                if value not in evidence:
                    evidence += [value]
                    counts += [0]
            counts = [x for _,x in sorted(zip(evidence,counts))]
            plt.bar(X_axis+(idx*0.2),counts,0.2,color=colors[idx],label=cat_lab[idx])
        scores.sort()
        plt.xticks(X_axis, scores, rotation=60)
        plt.yscale('log')
        plt.legend(prop={'size': 6})
        plt.ylabel('Number of LEEs')
        plt.xlabel('Total Score')
        plt.show()
        #plt.savefig('Output_Overview.png',bbox_inches = "tight",dpi=200)
        plt.close()

        # If Kind Score identifies sub-categories, add additional plot;
        #check with contradictions
        cat_vals = []
        for each in ['dir contradiction','sign contradiction','att contradiction']: cat_vals += [kind_values[each]]
        if len(set(cat_vals))>1:
            kind = list(kept['Kind Score'])

            strong_corr = kind.count(kind_values['strong corroboration'])
            weak_corr1 = kind.count(kind_values['weak corroboration1'])
            weak_corr2 = kind.count(kind_values['weak corroboration2'])
            weak_corr3 = kind.count(kind_values['weak corroboration3'])
            corrs = np.array([strong_corr, weak_corr1, weak_corr2, weak_corr3])
            mylabels = ["Strong Corroborations: "+str(strong_corr), "Weak Corroborations1: "+str(weak_corr1),
                        "Weak Corroborations2: "+str(weak_corr2), "Weak Corroborations3: "+str(weak_corr3)]
            mycolors = ['#235490','#2B65AD','#718EC5','#B0BEDA']
            plt.figure(figsize=(16, 4))
            plt.subplot(1, 4, 1)
            plt.pie(corrs,colors=mycolors)
            plt.legend(labels=mylabels,bbox_to_anchor=(0.2,0), loc="lower center",
                       bbox_transform=plt.gcf().transFigure)

            full_ext = kind.count(kind_values['full extension'])
            hang_ext = kind.count(kind_values['hanging extension'])
            int_ext = kind.count(kind_values['internal extension'])
            spec = kind.count(kind_values['specification'])
            exts = np.array([full_ext, hang_ext, int_ext, spec])
            mylabels = ["Full Extensions: "+str(full_ext), "Hanging Extensions: "+str(hang_ext),
                        "Internal Extensions: "+str(int_ext), "Specifications: "+str(spec)]
            mycolors = ['#24552D','#49A155','#7DBA84','#B6D5B8']
            plt.subplot(1, 4, 2)
            plt.pie(exts,colors=mycolors)
            plt.legend(labels=mylabels,bbox_to_anchor=(0.4,0), loc="lower center",
                       bbox_transform=plt.gcf().transFigure)

            dir_cont = kind.count(kind_values['dir contradiction'])
            sign_cont = kind.count(kind_values['sign contradiction'])
            att_cont =  kind.count(kind_values['att contradiction'])
            conts = np.array([dir_cont, sign_cont, att_cont])
            mylabels = ["Direction Contradictions: "+str(dir_cont), "Sign Contradictions: "+str(sign_cont),
                        "Attribute Contradictions: "+str(att_cont)]
            mycolors = ['#BE9735','#E2B441','#F9E6A9']
            plt.subplot(1, 4, 3)
            plt.pie(conts,colors=mycolors)
            plt.legend(labels=mylabels,bbox_to_anchor=(0.62,0), loc="lower center",
                       bbox_transform=plt.gcf().transFigure)

            flg1 = kind.count(kind_values['flagged1'])
            flg2 = kind.count(kind_values['flagged2'])
            flg3 = kind.count(kind_values['flagged3'])
            flgds = np.array([flg1, flg2, flg3])
            mylabels = ["Flagged 1: "+str(flg1), "Flagged 2: "+str(flg2),
                        "Flagged 3: "+str(flg3)]
            mycolors = ['#AA5626','#C7652E','#F1C6A8']
            plt.subplot(1, 4, 4)
            plt.pie(flgds,colors=mycolors)
            plt.legend(labels=mylabels,bbox_to_anchor=(0.82,0), loc="lower center",
                       bbox_transform=plt.gcf().transFigure)
            plt.savefig('Subcategory_Overview.png',bbox_inches = "tight",dpi=200)
            plt.close

    # If only visualizing one category:
    else:
        category = file_name.split('.')[0].split('_')[-1]
        mycolors = {"corroborations" : 'royalblue',
                    "extensions" : 'limegreen',
                    "contradictions" : 'gold',
                    "flagged" : 'darkorange'}

        plt.figure(figsize=(12, 6))
        #Evidence Score plots
        evidence = kept['Evidence Score'].value_counts().keys().tolist()
        counts = kept['Evidence Score'].value_counts().tolist()
        plt.subplot(2, 2, 1)
        plt.bar(evidence,counts,0.2,color=mycolors[category],label=category)
        plt.yscale('log')
        plt.legend(prop={'size': 6})
        plt.ylabel('Number of LEEs')
        plt.xlabel('Evidence Score')

        #Match Score plots
        plt.subplot(2, 2, 2)
        scores = [match_values['source present'],match_values['target present'],
                match_values['both present'],match_values['neither present']]
        scores.sort()
        X_axis = np.arange(len(scores))
        evidence = kept['Match Score'].value_counts().keys().tolist()
        counts = kept['Match Score'].value_counts().tolist()
        for value in scores:
            if value not in evidence:
                evidence += [value]
                counts += [0]
        counts = [x for _,x in sorted(zip(evidence,counts))]
        plt.bar(X_axis,counts,0.2,color=mycolors[category],label=category)
        plt.xticks(X_axis, scores)
        plt.legend(prop={'size': 6})
        plt.ylabel('Number of LEEs')
        plt.xlabel('Match Score')
        ticks_loc = X_axis


        #Total Score plots
        plt.subplot(2, 2, 3)

        scores = list(set(kept['Total Score']))
        scores = [float(i) for i in scores]
        X_axis = np.arange(len(scores))
        evidence = kept['Total Score'].value_counts().keys().tolist()
        counts = kept['Total Score'].value_counts().tolist()
        plt.bar(evidence,counts,0.2,color=mycolors[category],label=category)
        # plt.xticks(X_axis, evidence, rotation=45)
        plt.yscale('log')
        plt.legend(prop={'size': 6})
        plt.ylabel('Number of LEEs')
        plt.xlabel('Total Score')

        #Category Breakdown
        plt.subplot(2, 2, 4)
        #Pie Chart
        kind = list(kept['Kind Score'])
        cats = {'corroborations':['strong corroboration', 'weak corroboration1', 'weak corroboration2', 'weak corroboration3'],
                'corroborations colors':['#235490','#2B65AD','#718EC5','#B0BEDA'],
                'extensions':['full extension', 'hanging extension', 'internal extension', 'specification'],
                'extensions colors': ['#24552D','#49A155','#7DBA84','#B6D5B8'],
                'contradictions':['dir contradiction', 'sign contradiction', 'att contradiction','none'],
                'contradictions colors':['#BE9735','#E2B441','#F9E6A9','none'],
                'flagged':['flagged1', 'flagged2', 'flagged3', 'none'],
                'flagged colors':['#AA5626','#C7652E','#F1C6A8','none']}
        cat_df = pd.DataFrame(cats)
        #1: Get subcategory names
        sub_cats = cats[category]
        if 'none' in sub_cats: sub_cats.remove('none')
        mycolors = cats[category+' colors']
        if 'none' in mycolors: mycolors.remove('none')
        #2: see if separated (set of category values >2)
        cat_vals = []
        for each in sub_cats: cat_vals += [kind_values[each]]
        #3: If separated, get subcategory names
        if len(set(cat_vals))>1:
            mylabels = []
            counts = []
            for x in sub_cats:
                mylabels += [x+": "+str(kind.count(kind_values[x]))]
                counts += [kind.count(kind_values[x])]
            numbs = np.array(counts)
            plt.pie(numbs,colors=mycolors)
            plt.legend(labels=mylabels,bbox_to_anchor=(1,0.25), loc="lower right",
                    bbox_transform=plt.gcf().transFigure)
        else:
            mylabels = [category]
            counts = [output.shape[0]]
            numbs = np.array(counts)
            plt.pie(numbs,colors=[mycolors[1]],autopct=lambda p: '{:.0f}'.format(p * counts[0] / 100))
        plt.savefig(category+'_Overview.png',bbox_inches = "tight",dpi=200)
        plt.close
    return
