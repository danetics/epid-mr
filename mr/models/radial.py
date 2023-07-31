'''
Flag outliers based on Q statistics from radial regression
'''

import pandas as pd
import numpy as np
from scipy.stats import chi2

# Logic here is based on Bowden 2018 and the source code for the accompanying RadialMR R package (lines are permalinked)
# The RadialMR package uses iterative radial regressions, but this is not mentioned in the originating paper, so it is not implemented here 
# Many definitions of weights can be used, but using 'first order' inverse variance here throughout for simplicity
# Q statistics follow a chi2 distribution, and p-values are attained using R's pchisq(lower.tail=FALSE); inverse CDF
# chi2.sf() == survival function == 1-CDF is the equivalent: https://albertotb.com/Equivalence-between-distribution-functions-in-R-and-Python/ 

def flag_failures_cochran(data: pd.DataFrame, 
                          ivwrad_res: statsmodels.regression.linear_model.RegressionResults):
    '''
    Calculate Cochran's Q from radial IVW
    Return (i) a boolean mask which flags variant outliers, (ii) the summary Q statistic
    '''
    # Bowden 2018, Box 4 equation 2
    cochranq = (data['ratio_se']**-2) * (data['ratio']-ivwrad_res.params['1/ratio_se'])**2
    # https://github.com/WSpiller/RadialMR/blob/abe00170076284cb79ae711c96d3832da3879267/R/ivw_radial.R#L402
    cochranp = chi2.sf(cochranq, df=1) # one-sided
    return np.where((cochranp < 0.05), True, False), cochranq.sum()

def flag_failures_rucker(data: pd.DataFrame,
                         eggrad_res: statsmodels.regression.linear_model.RegressionResults):
    '''
    Calculate Rucker's Q from radial Egger
    Return (i) a boolean mask which flags variant outliers, (ii) the summary Q statistic
    '''
    # Bowden 2018, Box 4 equation 4
    ruckerq = (data['ratio_se']**-2) * (data['ratio']-(eggrad_res.params['Intercept']/data['ratio_se']**-1)-eggrad_res.params['1/ratio_se'])**2
    ruckerp = chi2.sf(ruckerq, df=1) # as for Cochran
    return np.where((ruckerp < 0.05), True, False), ruckerq.sum()

# chi2.sf(ruckerq.sum(), df=len(data)-1
    
    
    