'''
Flag outliers based on Q statistics from radial regression
'''

import statsmodels
import pandas as pd
import numpy as np
from scipy.stats import chi2

# Logic here is based on Bowden 2018 and the source code for the accompanying RadialMR R package (lines are permalinked)
# The RadialMR package uses iterative radial regressions, but this is not mentioned in the originating paper, so it is not implemented here 
# Many definitions of weights can be used, but using 'first order' inverse variance here throughout for simplicity
# Q statistics follow a chi2 distribution, and p-values are attained using R's pchisq(lower.tail=FALSE); inverse CDF
# chi2.sf() == survival function == 1-CDF is the equivalent: https://albertotb.com/Equivalence-between-distribution-functions-in-R-and-Python/ 

def calc_cochranq_per_variant(data: pd.DataFrame, fitted_beta: float, return_pvals=True):
    '''
    Calculate row-wise Cochran's Q relative to a fitted regression slope
    Can be used for radial and normal Q-statistics: all that changes is model beta (summary effect size term)
    '''
    # Bowden 2018, Box 4 equation 2: normal Q statistic (weighted sum of squared mean differences)
    # https://github.com/WSpiller/RadialMR/blob/abe00170076284cb79ae711c96d3832da3879267/R/ivw_radial.R#L402
    cochranq = (data['ratio_se']**-2) * (data['ratio']-fitted_beta)**2
    cochranp = chi2.sf(cochranq, 1)
    if return_pvals:
        return cochranq, cochranp
    return cochranq

def calc_ruckerq_per_variant(data: pd.DataFrame, eggrad_res):
    '''
    Calculate row-wise Rucker's Q from radial Egger
    '''
    # Bowden 2018, Box 4 equation 4
    ruckerq = (data['ratio_se']**-2) * (data['ratio']-(eggrad_res.params['Intercept']/data['ratio_se']**-1)-eggrad_res.params['ratio_inv_se'])**2
    ruckerp = chi2.sf(ruckerq, 1)
    return ruckerq, ruckerp

def calc_summary_heterogeneity_statistics(qstats: pd.Series):
    '''
    Given an array of variant-wise Q statistics, return the summary Q statistics, its p-value and corresponding I2
    '''
    q = qstats.sum()
    pval = chi2.sf(q, len(qstats)-1) # one-sided
    i2 = 100*((q-(len(qstats)-1))/q)
    return q, pval, i2

    
    
    