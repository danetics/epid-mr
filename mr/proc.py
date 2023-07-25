'''
Subfunctions to prepare top-level input data for univariable and multivariable analysis wrappers
'''

# TODO - aligned freq
# TODO - I2 statistic 
# TODO - proxy expansion 

import pandas as pd
import numpy as np
    
def align_alleles(data: pd.DataFrame) -> pd.DataFrame:
    ''' 
    Ensure that effect alleles match between the studies 
    '''
    data = data[(((data['effect_allele_x']==data['effect_allele_y']) & 
                  (data['other_allele_x']==data['other_allele_y'])) | 
                 ((data['effect_allele_x']==data['other_allele_y']) & 
                  (data['other_allele_x']==data['effect_allele_y'])))]
    data['flip_beta_y'] = np.where(((data['effect_allele_x']==data['other_allele_y']) & 
                                  (data['other_allele_x']==data['effect_allele_y'])), True, False)
    data['beta_y'] = np.where(data['flip_beta_y'], -(data['beta_y']), data['beta_y'])
    return data

def align_betas(data: pd.DataFrame) -> pd.DataFrame:
    ''' 
    Align outcome betas to exposure-increasing allele
    '''
    # Betas are already aligned to the same effect allele, so we just need to flip them as required
    data['beta_y'] = np.where(data['beta_x']>=0, data['beta_y'], -(data['beta_y']))
    data['beta_x'] = np.where(data['beta_x']>=0, data['beta_x'], -(data['beta_x']))
    data['ratio'] = data['beta_y']/data['beta_x']
    return data

def keep_reqcols(data):
    return data[['chr','bp','beta_x','se_x','pval_x','n_x','beta_y','se_y','pval_y','n_y']]

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
    out = keep_reqcols(align_betas(align_alleles(data)))
    out['chrpos'] = build_chrpos(data['chr'], data['bp'])
    out['ratio'] = calc_ratio(out['beta_x'], out['beta_y'])
    out['weight'] = calc_inverse_se2(out['se_y'])
    return out 




    
    

    
    