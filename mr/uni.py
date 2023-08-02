'''
Top-level module to run univariaable MR models
'''
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from scipy.stats import chi2
from .models import med, steig, het
from .plot import panel
from sumstats.extract import instrument

def run_main(data: pd.DataFrame) -> dict:
    return {
        'ivw': smf.wls('beta_y ~ beta_x - 1', data, weights=data['se_y']**-2).fit(),
        'egger': smf.wls('beta_y ~ beta_x', data, weights=data['se_y']**-2).fit(),
        'wm': med.run_wm(data)
    }

def run_steiger(data: pd.DataFrame): 
    data['steiger_fail'] = steig.flag_failures(data)
    steiger_filt = data[~(data['steiger_fail'])]
    res = {
        'ivw_steig_filtered': smf.wls('beta_y ~ beta_x - 1', steiger_filt, weights=steiger_filt['se_y']**-2).fit()
    }
    return data, res

def run_radial(data: pd.DataFrame):
    res = {
        'radial': smf.ols('ratio_z ~ ratio_inv_se - 1', data).fit(),
        'egger_radial': smf.ols('ratio_z ~ ratio_inv_se', data).fit()
    }
    data['cochranq_radial'], cochradp = het.calc_cochranq_per_variant(data, res['radial'].params['ratio_inv_se'])
    data['ruckerq_radial'], ruckradp = het.calc_ruckerq_per_variant(data, res['egger_radial'])
    data['radial_fail'] = np.where(((cochradp < 0.05) | (ruckradp < 0.05)), True, False)
    radial_filt = data[~data['radial_fail']]
    filt = {
        'ivw_radial_filtered': smf.wls('beta_y ~ beta_x - 1', radial_filt, weights=radial_filt['se_y']**-2).fit(),
        'egger_radial_filtered': smf.wls('beta_y ~ beta_x', radial_filt, weights=radial_filt['se_y']**-2).fit()
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
    data['cochranq_ivw'] = het.calc_cochranq_per_variant(data, main['ivw'].params['beta_x'], return_pvals=False)
    data['cochranq_ivw_radial_filtered'] = het.calc_cochranq_per_variant(data[~data['radial_fail']], radial['ivw_radial_filtered'].params['beta_x'], return_pvals=False)
    # Rucker Q could also be calculated here for main and radial filtered models, but not bothering
    return data, {**main, **steiger, **radial}

def from_hdfpaths(xpath: str, 
                  ypath: str, 
                  get_proxies=True,
                  omitx = None,
                  rsidsx = None,
                  omity = None,
                  rsidsy = None):
    # When there are multiple signals or instrument tables we can also add in an option for that, to pass to sumstats.extract.instrument
    data_xy = instrument.get_analytic_dataframe(xpath, ypath, get_proxies, omitx, rsidsx)
    data_yx = instrument.get_analytic_dataframe(ypath, xpath, get_proxies, omity, rsidsy)
    proc_xy, res_xy = run_analyses(data_xy)
    proc_yx, res_yx = run_analyses(data_yx)
    panel.mr_panel(proc_xy, proc_yx, res_xy, res_yx)
    
    

    
    
        
    