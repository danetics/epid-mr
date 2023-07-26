'''
Subfunctions for calculating and applying MR-Steiger based identification of outlying variants
'''
import pandas as pd
import numpy as np
from scipy.stats import chi2
from . import reg

# NOTE - other definitions of first-order weights can be used, but using inverse variance here for simplicity and consistency
# NOTE - Q statistics are treated as Chi2
# NOTE - the ratio hereis synonymous with Bj in the paper - the effect size for the variant

def calc_cochranq(ratio: pd.Series, 
                  weight: pd.Series) -> pd.Series:
    ''' Run MR-IVW on radial scale, and use to calculate per-variant Cochran's Q '''
    # sqrt(1/se**2) == 1/se
    # Below is (ratio/se) ~ 1/se without an intercept, per Bowden 2018, Box 4, equation 1
    radialreg = reg.constrained_wls(xvar = np.sqrt(weight),
                                    yvar = (ratio * np.sqrt(weight)))
    # Bowden 2018, Box 4 equation 2
    return weight*(ratio-radialreg['beta'])**2

def calc_ruckerq(ratio: pd.Series, 
                 weight: pd.Series) -> pd.Series:
    ''' Run MR-Egger on radial scale, and use to calculate per-variant Rucker's Q '''
    # As for Cochran but with an intercept: Bowden 2018 Box 4 equation 3
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
    cochranp = chi2.sf(cochranq, df=1) # one-sided
    outliers = data[cochranp <= 0.05]
    summary = {
        'dropped': outliers['chrpos'].tolist(),
        'dropped_n': len(outliers),
        'qstat': round(cochranq.sum(), 2),
        'qstat_pval': chi2.sf(cochranq.sum(), df=len(data)-1)
    }
    return (data[cochranp > 0.05], summary)

def apply_radial_egger(data: pd.DataFrame) -> tuple:
    '''
    '''
    ruckerq = calc_ruckerq(data['ratio'], data['weight'])
    ruckerp = chi2.sf(ruckerq, df=1) # one-sided
    outliers = data[ruckerp <= 0.05]
    summary = {
        'dropped': outliers['chrpos'].tolist(),
        'dropped_n': len(outliers),
        'qstat': round(ruckerq.sum(), 2),
        'qstat_pval': chi2.sf(ruckerq.sum(), df=len(data)-1)
    }
    return (data[ruckerp > 0.05], summary)