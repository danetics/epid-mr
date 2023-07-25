'''
Top-level wrappers to run all univariable MR methods and sensitivity analyses
'''

import pandas as pd
import hdf, proc, reg, med, rad, steig

def format_output(data, xpath, ypath) -> pd.DataFrame:
    data.insert(1, 'outcome', hdf.read_metadata(ypath)['trait'])
    data.insert(1, 'exposure', hdf.read_metadata(xpath)['trait'])
    data = data.rename(columns={'index': 'method'})
    return data

def run_main(data: pd.DataFrame) -> pd.DataFrame:
    ''' Run univariable MR analyses using all methods, without filtering '''
    res = {
        'MR-IVW': reg.constrained_wls(data['beta_x'], data['beta_y'], weights=data['weight']),
        'MR-Egger': reg.unconstrained_wls(data['beta_x'], data['beta_y'], weights=data['weight']),
        'MR-WM':med.run_wm(data['ratio'], data['weight'])
    }
    return pd.DataFrame.from_dict(res, orient='index').reset_index()
                 
def run_steiger(data):
    ''' Run regression-based MRs using only variants passing MR-Steiger filtering '''
    steiger = steig.apply_steiger(data)
    res = {
        'MR-IVW': (reg.constrained_wls(steiger[0]['beta_x'], steiger[0]['beta_y'], weights=(steiger[0]['weight']))| steiger[1]),
        'MR-Egger': (reg.unconstrained_wls(steiger[0]['beta_x'], steiger[0]['beta_y'], weights=(steiger[0]['weight'])) | steiger[1])
    }
    return pd.DataFrame.from_dict(res, orient='index').reset_index()

def run_radial(data):
    ''' Run regression-based MRs using only variants passing MR-Radial filtering'''
    radivw = rad.apply_radial_ivw(data)
    radegg = rad.apply_radial_egger(data)
    res = {
        'MR-IVW': (reg.constrained_wls(radivw[0]['beta_x'], radivw[0]['beta_y'], weights=(radivw[0]['weight'])) | radivw[1]),
        'MR-Egger': (reg.unconstrained_wls(radegg[0]['beta_x'], radegg[0]['beta_y'], weights=(radegg[0]['weight'])) | radegg[1])
    }
    return pd.DataFrame.from_dict(res, orient='index').reset_index()

def run_analyses(xpath: str, ypath: str) -> pd.DataFrame:
    '''
    '''
    instrument = hdf.load_signals(xpath)
    data = hdf.extract_instrument(ypath, instrument)
    data = proc.preprocess_data(data)
    return {
        'main': format_output(run_main(data), xpath, ypath), 
        'steiger': format_output(run_steiger(data), xpath, ypath), 
        'radial': format_output(run_radial(data), xpath, ypath)
    }
        

#----->>>>>----->>>>> MAIN <<<<<-----<<<<<-----

def run_mr(xpath: str,
           ypath: str,
           bidirectional=False):
    '''
    '''
    res_xy = run_analyses(xpath, ypath)
    if bidirectional:
        res_yx = run_analyses(ypath, xpath)
        return (res_xy, res_yx) 
    return res_xy
               
    






               