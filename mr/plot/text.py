'''
Generate text to be included in plot area
'''

import pandas as pd
from ..models import het

def get_variant_counts(data: pd.DataFrame, highlight: str) -> str:
    '''
    Given an analytic dataframe, return the counts of variants therein as a formatted string for display on the plot
    If a highlight pd.Series is also passed, remove those variants from the data before counting
    '''
    # Use case for removing the highlighted variants is the filtered models
    filtered = data[~data[highlight]]
    return f"n={len(data)} ({len(data[data['rsid_x']!=data['rsid_y']])} proxies) subset={len(filtered)} ({round((len(filtered)/len(data))*100, 0)}%)"

def get_qstats(data, reglines: list):
    if 'ivw' in reglines:
        q, pval, i2 = het.calc_summary_heterogeneity_statistics(data['cochranq_ivw'])
    elif 'ivw_radial_filtered' in reglines:
        q, pval, i2 = het.calc_summary_heterogeneity_statistics(data['cochranq_ivw_radial_filtered'])
    elif 'radial' in reglines:
        q, pval, i2 = het.calc_summary_heterogeneity_statistics(data['cochranq_radial'])
    cochran = f"Cochran's Q p={pval:.1e}, (I\u00B2={round(i2,1)})"
    if 'egger_radial' in reglines:
        q, pval, i2 = het.calc_summary_heterogeneity_statistics(data['ruckerq_radial'])
        rucker = f"Rucker's Q p={pval:.1e}, (I\u00B2={round(i2,1)})"
        return f"{cochran}\n{rucker}"
    else:
        return cochran