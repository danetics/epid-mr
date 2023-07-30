'''
Top-level module to run univariaable MR models
'''

import statsmodels.formula.api as smf


def run_main(data: pd.DataFrame) -> dict:
    return {
        'ivw': smf.wls('beta_y ~ beta_x - 1', data, weights=data['se_y']**-2).fit(),
        'egger': smf.wls('beta_y ~ beta_x', data, weights=data['se_y']**-2).fit(),
        'wm': models.med.run_wm(data['ratio'], weights=data['ratio_se']**-2)
    }

def run_steiger(data: pd.DataFrame): 
    data['steiger_fail'] = models.steig.flag_failures(data)
    if any(data['steiger_fail']):
        if len(data[~data['steiger_fail']]) >= 5:
            res = {
                'ivw_steig_filtered': smf.wls('beta_y ~ beta_x - 1', data, subset=(~data['steiger_fail']), weights=data['se_y']**-2).fit(),
            }
            return data, res
    return data, dict()

def run_radial(data: pd.DataFrame):
    data['radial_fail'] = models.radial.flag_failures(data)
    res = {
        'radial': smf.ols('ratio/ratio_se ~ (1/ratio_se) - 1', data, subset=(~data['radial_fail'])).fit(),
        'egger_radial': smf.ols('ratio/ratio_se ~ (1/ratio_se)', data, subset=(~data['radial_fail'])).fit()
    }
    if any(data['radial_fail']):
        if len(data[~data['radial_fail']]) >= 5:
            filt = {
                'ivw_radial_filtered': smf.wls('beta_y ~ beta_x - 1', data, subset=(~data['radial_fail']), weights=data['se_y']**-2).fit()
                'egger_radial_filtered': smf.wls('beta_y ~ beta_x', data, subset=(~data['radial_fail']), weights=data['se_y']**-2).fit()
            }
            return data, {**res, **filt}
    return data, res

def run_analyses(data):
    '''
    For a given preprocessed dataframe, add Steiger and Radial flags and run all models
    Return the modified data and a dict of results
    '''
    main = run_main(data)
    data, steiger = run_steiger(data)
    data, radial = run_radial(data)
    return data, {**main, **steiger, **radial}

    
    
        
    