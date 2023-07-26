'''
Subfunctions to prepare top-level input data for univariable and multivariable analysis wrappers
'''

# TODO - I2 statistic

import pandas as pd
import numpy as np



def keep_reqcols(data):
    return data[['rsid_x','rsid_y','beta_x','se_x','pval_x','n_x','entity_x','beta_y','se_y','pval_y','n_y','entity_y']]

def build_chrpos(chr: pd.Series,
                 bp: pd.Series) -> pd.Series:
    ''' Add a composite marker string for easy list of variants etc '''
    return pd.Series('chr' + chr.astype(str) + ':' + bp.astype(str))

def calc_inverse_se2(se_y: pd.Series):
    return se_y**-2

def calc_ratio(beta_x: pd.Series, beta_y: pd.Series) -> pd.Series:
    return beta_y/beta_x

def preprocess_data(data: pd.DataFrame) -> pd.DataFrame:
    ''' Wrapper for data processing steps '''
    out = keep_reqcols(data)
    assert not any(out.duplicated('rsid_x'))
    assert not any(out.duplicated('rsid_y')) 
    assert not any(out.duplicated(['rsid_x','rsid_y']))
    out['ratio'] = calc_ratio(out['beta_x'], out['beta_y'])
    out['weight'] = calc_inverse_se2(out['se_y'])
    return out 




    
    

    
    