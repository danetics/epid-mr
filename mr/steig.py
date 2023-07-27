'''
Functions to calculate Steiger-MR Z scores and define Steiger failures
Logic is based on original Hemani 2017 paper, plus associated R package
'''

import pandas as pd
import numpy as np
from scipy.stats import norm, f

# NOTE: Steiger calculations are dependent on per variant pval aone, which makes them sensitive to pval type precision

def calc_f_from_sumstats(pval: pd.Series, 
                         n: pd.Series) -> pd.Series:
    ''' Given the pvalue and n of a SNP in sumstats, calculate its F-statistic '''
    # Use the F-distribution inverse survival function
    return f.isf(pval, 1, n-1)

def calc_r_from_f(f: pd.Series, 
                  n: pd.Series) -> pd.Series:
    ''' Convert F to correlation r (unsquared) '''
    return np.sqrt(f/(n-2+f))

def transform_r_to_zg(r: pd.Series) -> pd.Series:
    ''' Use Fisher transformation to convert per-variant r to per-variant Zg '''
    # Hemani 2017 - p16, equation 2 (Fisher's Z-transformation)
    return 0.5 * np.log((1+r)/(1-r)) # np.log is ln() 

def calc_steiger_z(zgx: pd.Series, 
                   zgy: pd.Series,
                   nx: pd.Series,
                   ny: pd.Series) -> pd.Series:
    ''' Calculate Steiger's Z using formula for 2SMR '''
    # Hemani 2017 - p17, equation 1 (Steiger Z in two-sample MR)
    return (zgx-zgy)/np.sqrt((1/(nx-3)) + (1/(ny-3)))

def apply_steiger(data: pd.DataFrame) -> tuple:
    '''
    '''
    zgx = transform_r_to_zg(calc_r_from_f(calc_f_from_sumstats(data['pval_x'], data['n_x']), data['n_x']))
    zgy = transform_r_to_zg(calc_r_from_f(calc_f_from_sumstats(data['pval_y'], data['n_y']), data['n_y']))
    steigerz = calc_steiger_z(zgx, zgy, data['n_x'], data['n_y'])
    steigerp = norm.sf(steigerz)
    failures = data[((steigerp <= 0.05) & (steigerz < 0))]
    summary = {
        'dropped': failures['rsid_x'].tolist(),
        'dropped_n': len(failures)
    }
    return (data[~data.index.isin(failures.index)], summary)