'''
Subfunctions for calculating and applying MR-Steiger based identification of outlying variants
'''

import pandas as pd
import numpy as np
from scipy.stats import chi2
from . import reg

# Logic here is based on Bowden 2018 and the source code for the accompanying RadialMR R package (lines are permalinked)
# The RadialMR package uses iterative radial regressions, but this is not mentioned in the originating paper, so it is not implemented here 
# Many definitions of weights can be used, but using 'first order' inverse variance here throughout for simplicity
# Q statistics follow a chi2 distribution, and p-values are attained using R's pchisq(lower.tail=FALSE); inverse CDF
# chi2.sf() == survival function == 1-CDF is the equivalent, see https://albertotb.com/Equivalence-between-distribution-functions-in-R-and-Python/ 
# the ratio measure used here == Bj in the paper == per-variant effect size

def calc_cochranq(ratio: pd.Series, 
                  weight: pd.Series) -> pd.Series:
    ''' 
    Run MR-IVW on radial scale, and use to calculate per-variant Cochran's Q 
    '''
    # Precalculated weights are 1/se**2.  Paper uses 1/se to weight.  sqrt(1/se**2) == 1/se
    # Below regression is ((ratio/se) ~ 1/se - 1), from Bowden 2018, Box 4, equation 1
    radialreg = reg.constrained_wls(xvar = np.sqrt(weight),
                                    yvar = (ratio * np.sqrt(weight)))
    # Bowden 2018, Box 4 equation 2
    return weight*(ratio-radialreg['beta'])**2

def calc_ruckerq(ratio: pd.Series, 
                 weight: pd.Series) -> pd.Series:
    ''' Run MR-Egger on radial scale, and use to calculate per-variant Rucker's Q '''
    # As for Cochran but with intercept: Bowden 2018 Box 4 equation 3
    radialreg = reg.unconstrained_wls(xvar = np.sqrt(weight),
                                      yvar = (ratio * np.sqrt(weight)))
    # Bowden 2018 Box 4 equation 4
    return weight*(ratio-(radialreg['intercept']/np.sqrt(weight))-radialreg['beta'])**2

def apply_radial_ivw(data: pd.DataFrame) -> tuple:
    '''
    Calculate per-variant Q-stat, infer pval and use to define radial outliers
    Return subset of non-outlier variants and a dict of summary information
    '''
    cochranq = calc_cochranq(data['ratio'], data['weight'])
    # Individual contribution to heterogeneity measured by pval of Qj
    # https://github.com/WSpiller/RadialMR/blob/abe00170076284cb79ae711c96d3832da3879267/R/ivw_radial.R#L402
    cochranp = chi2.sf(cochranq, df=1) # one-sided
    outliers = data[cochranp < 0.05]
    summary = {
        'dropped': outliers['rsid_x'].tolist(),
        'dropped_n': len(outliers),
        'qstat': cochranq.sum(),
        'qstat_pval': chi2.sf(cochranq.sum(), df=len(data)-1) # one-sided
    }
    return (data[~(data['rsid_x'].isin(outliers['rsid_x']))], summary)

def apply_radial_egger(data: pd.DataFrame) -> tuple:
    '''
    '''
    ruckerq = calc_ruckerq(data['ratio'], data['weight'])
    # Get pval as for Cochran
    ruckerp = chi2.sf(ruckerq, df=1) # one-sided
    outliers = data[ruckerp < 0.05]
    summary = {
        'dropped': outliers['rsid_x'].tolist(),
        'dropped_n': len(outliers),
        'qstat': ruckerq.sum(),
        'qstat_pval': chi2.sf(ruckerq.sum(), df=len(data)-1) # one-sided
    }
    return (data[~(data['rsid_x'].isin(outliers['rsid_x']))], summary)