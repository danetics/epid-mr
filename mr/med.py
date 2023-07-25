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

def bootstrapped_se(ratio: pd.Series,
                    weight: pd.Series):
    # Can simply bootstrap a vector of (ratio * weight) and calculate SE        
    return bootstrap(((ratio*weight),), np.median, method='percentile')

# ----->>>>> Function 

def run_wm(ratio: pd.Series,
           weight: pd.Series) -> dict:
    ''' 
    Weighted median MR with specified weights
    '''
    out = {
        'beta': robustats.weighted_median(ratio, weight),
        'se': bootstrapped_se(ratio, weight).standard_error,
        'nvar': len(ratio)
    }
    out['pval'] = 2 * norm.sf(np.abs(out['beta']/out['se'])) 
    return out

