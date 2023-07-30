'''
Subfunctions to calculate median (non-regression) based MR estimates
Functions to penalise first-order weights have been removed for now given limited utility of PWM
'''

import pandas as pd
import numpy as np
import robustats
from scipy.stats import norm, bootstrap
from statsmodels.formula import api as smf

# ----->>>>> Subfunctions

def get_weighted_median(sorted: pd.DataFrame) -> float:
    # Easier to use full DF so that values are sorted together without additional processing
    sorted['cumulative_weight'] = np.cumsum(sorted['ratio_se']**-2)
    totalweight = sorted['cumulative_weight'].sum
    median_idx = np.searchsorted(sorted['ratio'], totalweight/2)
    if totalweight%2==1:
        return sorted['ratio'][median_idx]
    return (sorted['ratio'][median_idx-1] + sorted['ratio'][median_idx+1])/2
        
def bootstrapped_se(sorted: pd.DataFrame) -> float:
    # Can simply bootstrap the sorted vector of (ratio * weight) and calculate SE        
    res = bootstrap(((sorted['ratio'] * sorted['ratio_se']**-2),), np.median, method='percentile')
    return res.standarderror

# ----->>>>> Function 

def run_wm(data: pd.DataFrame) -> dict:
    ''' 
    Weighted median MR with specified weights
    '''
    sorted = data.sort_values('ratio')
    return {
        'beta': get_weighted_median(sorted),
        'se': bootstrapped_se(sorted),
        'nvar': len(data),
        'pval': 2 * norm.sf(np.abs(get_weighted_median(sorted)/bootstrapped_se(sorted)))
    }
  

