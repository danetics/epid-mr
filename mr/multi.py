'''
Top-level wrapper to run all multivariable MR methods and sensitivity analyses
'''

import pandas as pd
import hdf, proc, reg, rad, steig

def run_ivw(x1, x2, y, weights):
    '''
    Multivariable MR-IVW with one covariate, implemented as a two-step regression
    '''
    stepone = reg.constrained_wls(xvar=x2, yvar=y, weights=weights)
    if not stepone=={}:
        residuals = stepone['resid']
        steptwo = reg.constrained_wls(xvar=x1, yvar=residuals)
    