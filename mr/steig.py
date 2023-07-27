'''
Functions to calculate Steiger-MR Z scores and define Steiger failures
'''

import pandas as pd
import numpy as np
from scipy.stats import norm, f

# Steiger calculations are dependent on per variant pval aone, which makes them sensitive to pval type precision
# Logic here is based on original Hemani 2017 paper and the worked example from that publication (code lines are permalinked)
# Equivalence between R and scipy.stats distribution function from https://albertotb.com/Equivalence-between-distribution-functions-in-R-and-Python/

def calc_f_from_sumstats(pval: pd.Series, 
                         n: pd.Series) -> pd.Series:
    ''' 
    Calculate F-statistic of a variant from pval and n in sumstats using the f distribution inverse survival function
    '''
    # https://github.com/explodecomputer/causal-directions/blob/916f528cbf8d1ccd9bc7eac9f33f6b658e683ded/scripts/mr_directionality_analysis.R#L39
    # f.isf == R qf(lower.tail=FALSE)
    return f.isf(pval, 1, n-1) # one sided 

def calc_r_from_f(f: pd.Series, 
                  n: pd.Series) -> pd.Series:
    ''' 
    Convert F to correlation r (unsquared) using accepted approach with 1 DF for the numerator
    '''
    # https://github.com/explodecomputer/causal-directions/blob/916f528cbf8d1ccd9bc7eac9f33f6b658e683ded/scripts/mr_directionality_analysis.R#L40
    # http://davidileitman.com/wp-content/uploads/2014/04/EffectSizeFormulas.pdf
    return np.sqrt(f/(f+(n-2)))

def transform_r_to_zg(r: pd.Series) -> pd.Series:
    ''' 
    Use Fisher transformation to convert per-variant r to per-variant Zg
    '''
    # Hemani 2017, p16, equation 2
    return 0.5 * np.log((1+r)/(1-r)) # np.log is ln() 

def calc_steiger_z(zgx: pd.Series, 
                   zgy: pd.Series,
                   nx: pd.Series,
                   ny: pd.Series) -> pd.Series:
    ''' 
    Calculate per-variant Steiger's Z using formula for 2SMR 
    '''
    # Hemani 2017, p17, equation 1
    return (zgx-zgy)/np.sqrt((1/(nx-3)) + (1/(ny-3)))

def calc_ratio_pval(sey, bx, ratio):
    '''
    Calculate wald ratio p-value for the by/bx ratio
    This is the per-variant p-value for a single-variant MR using that variant as the instrument
    '''
    # https://github.com/explodecomputer/causal-directions/blob/916f528cbf8d1ccd9bc7eac9f33f6b658e683ded/scripts/mr_directionality_analysis.R#L143
    # norm.sf == R pnorm(lower.tail=FALSE)
    wald_se = sey/abs(bx)
    wald_z = abs(ratio)/wald_se # ratio already pre-calculated in input
    return (2 * norm.sf(wald_z)) # two-sided
    
def get_failures(data, steigerz, waldp):
    '''
    Subset per-variant failures, where we accept y->x causation based on sign(steigerz) and the steigerz p-value
    '''
    # Hemani 2017, p16
    steigerp = 2 * norm.sf(steigerz)
    return data[((steigerp < 0.05) & (steigerz < 0) & (waldp < 0.05))]

def apply_steiger(data: pd.DataFrame) -> tuple:
    '''
    '''
    print(data.head())
    zgx = transform_r_to_zg(calc_r_from_f(calc_f_from_sumstats(data['pval_x'], data['n_x']), data['n_x']))
    zgy = transform_r_to_zg(calc_r_from_f(calc_f_from_sumstats(data['pval_y'], data['n_y']), data['n_y']))
    steigerz = calc_steiger_z(zgx, zgy, data['n_x'], data['n_y'])
    waldp = calc_ratio_pval(data['se_y'], data['beta_x'], data['ratio'])
    failures = get_failures(data, steigerz, waldp)
    summary = {
        'dropped': failures['rsid_x'].tolist(),
        'dropped_n': len(failures)
    }
    return (data[~(data['rsid_x'].isin(failures['rsid_x']))], summary)