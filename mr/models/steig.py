'''
Functions to calculate Steiger Z scores and define Steiger failures
'''

import pandas as pd
import numpy as np
from scipy.stats import norm, f

# Steiger calculations are dependent on per variant pval aone, which makes them sensitive to pval type precision
# Logic here is based on original Hemani 2017 paper and the worked example from that publication (code lines are permalinked)
# Equivalence between R and scipy.stats distributions from https://albertotb.com/Equivalence-between-distribution-functions-in-R-and-Python/

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

def flag_failures(data: pd.DataFrame) -> pd.Series:
    '''
    Identify variants where we accept a y -> x causal direction based on the Steiger logic
    Return a boolean mask which flags these Steiger failures
    '''
    # Calculate per-variant Steiger's Z
    zgx = transform_r_to_zg(calc_r_from_f(calc_f_from_sumstats(data['pval_x'], data['n_x']), data['n_x']))
    zgy = transform_r_to_zg(calc_r_from_f(calc_f_from_sumstats(data['pval_y'], data['n_y']), data['n_y']))
    steigerz = calc_steiger_z(zgx, zgy, data['n_x'], data['n_y'])
    
    # Define failures based on Hemani 2017, p16
    # Use steigerp, sign(steigerz) and the original Wald ratio pvalue
    steigerp = 2 * norm.sf(steigerz)
    ratiop = 2 * norm.sf(np.abs(data['ratio']/data['ratio_se']))
    return(np.where(((steigerp < 0.05) & (steigerz < 0) & (ratiop < 0.05)), True, False))
    
    