'''
Top-level module to run univariaable MR models
'''
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

def run_main(data: pd.DataFrame) -> dict:
    return {
        'ivw': smf.wls('beta_y ~ beta_x - 1', data, weights=data['se_y']**-2).fit(),
        'egger': smf.wls('beta_y ~ beta_x', data, weights=data['se_y']**-2).fit(),
        'wm': models.med.run_wm(data)
    }

def run_steiger(data: pd.DataFrame): 
    data['steiger_fail'] = models.steig.flag_failures(data)
    res = {
        'ivw_steig_filtered': smf.wls('beta_y ~ beta_x - 1', data, subset=(~data['steiger_fail']), weights=data['se_y']**-2).fit()
    }
    return data, res

def run_radial(data: pd.DataFrame):
    res = {
        'radial': smf.ols('ratio/ratio_se ~ (1/ratio_se) - 1', data).fit(),
        'egger_radial': smf.ols('ratio/ratio_se ~ (1/ratio_se)', data).fit()
    }
    data['radial_fail_coch'] = models.radial.flag_failures_cochran(data, res['radial'])
    data['radial_fail_ruck'] = models.radial.flag_failures_cochran(data, res['egger_radial'])
    filt = {
        'ivw_radial_filtered': smf.wls('beta_y ~ beta_x - 1', data, subset=(~data['radial_fail_coch']), weights=data['se_y']**-2).fit()
        'egger_radial_filtered': smf.wls('beta_y ~ beta_x', data, subset=(~data['radial_fail_ruck']), weights=data['se_y']**-2).fit()
    }
    return data, {**res, **filt}

def run_analyses(data):
    '''
    For a given preprocessed dataframe, add Steiger and Radial flags and run all models
    Return the modified data and a dict of results
    '''
    main = run_main(data)
    data, steiger = run_steiger(data)
    data, radial = run_radial(data)
    return data, {**main, **steiger, **radial}

    
    
        
    