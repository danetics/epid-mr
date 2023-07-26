'''
Subfunctions to run (weighted) OLS regressions
Either as a direct analysis (IVW, Egger) or as part of another process
'''

import pandas as pd
import statsmodels.formula.api as smf

def parse_results(res=None) -> dict:
    ''' 
    Parse or calculate required quantities from regression results
    Allows for one term and an optional intercept, returning a dict 
    '''
    # If the regression has not been run, return an empty dict 
    xname = res.params.keys()[0]
    term = {
        'beta': res.params[xname],
        'se': res.bse[xname],
        'pval': res.pvalues[xname],
        'nvar': int(res.nobs),
        'resid': res.resid
    }
    if 'Intercept' in res.params:
        intercept = {
            'intercept': res.params['Intercept'],
            'se_intercept': res.bse['Intercept'],
            'pval_intercept': res.pvalues['Intercept'],
        }
        return (term | intercept)
    return term

def prepare_olsdata(xvar: pd.Series,
                    yvar: pd.Series) -> pd.DataFrame:
    ''' Concatenate xvar, yvar series into analytic dataframe for OLS '''
    try:
        data = pd.concat([xvar, yvar], axis=1)
    except valueError:
        print('There are no intersecting datapoints for OLS: skipping')
        return pd.DataFrame()
    if not len(data)>5:
        print('There are <5 intersecting datapoints for OLS: skipping')
        return pd.DataFrame()
    return data

def constrained_wls(xvar: pd.Series,
                    yvar: pd.Series,
                    weights=None) -> dict:
    ''' 
    Run WLS regression with no intercept, return a dict of results
    Used directly for MR-IVW estimates
    '''
    data = prepare_olsdata(xvar, yvar)
    if data.empty: 
        return {}
    weights = pd.Series(1, index=range(len(xvar))) if weights is None else weights
    return parse_results(smf.wls('yvar ~ xvar - 1', data=data, weights=weights).fit())

def unconstrained_wls(xvar: pd.Series, 
                      yvar: pd.Series, 
                      weights=None) -> dict:
    ''' 
    Run WLS regression with an intercept, return a dict of results
    Used directly for MR-Egger estimates
    '''
    data = prepare_olsdata(xvar, yvar)
    if data.empty:
        return {}
    weights = pd.Series(1, index=range(len(xvar))) if weights is None else weights
    return parse_results(smf.wls('yvar ~ xvar', data=data, weights=weights).fit())

