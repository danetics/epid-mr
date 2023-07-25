'''
Top-level wrapper to run all multivariable MR methods and sensitivity analyses
'''

def run_main(data: pd.DataFrame) -> pd.DataFrame:
    ''' Run multivariable MR analyses using all methods, without filtering '''
    res = {
        'MR-IVW': reg.constrained_wls(data['beta_x'], 
                                      data['beta_y'],
                                      covar=data['beta_covar']
                                      weights=data['weight']),
        'MR-Egger': reg.unconstrained_wls(data['beta_x'], data['beta_y'], weights=data['weight'])
    }
    return pd.DataFrame.from_dict(res, orient='index').reset_index()