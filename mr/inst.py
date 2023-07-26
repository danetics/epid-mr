import pandas as pd
import numpy as np

def extract_variants(variants: list, store: pd.HDFStore) -> pd.DataFrame:
    ''' 
    For a list of variants as tuples [(chr, rsid)],
    Extract variants from sumstats/ group of the provided HDFStore 
    '''
    extracted = list()
    for x in variants:
        data = store.select(f'/sumstats/chr{x[0]}', where='rsid=x[1]')
        extracted = extracted + [data]
    return pd.concat(extracted)

def find_proxies(variants: list, proxymap: pd.DataFrame) -> pd.DataFrame:
    '''
    For a list of variants as tuples[(chr, bp_b37, rsid)]
    Extract possible proxies from the proxies/main table of the HDFStore
    '''
    vardf = pd.DataFrame(variants, columns=['chr','rsid'])
    return vardf.merge(proxymap, left_on='rsid', right_on='rsid_a', how='inner')

def keep_best_proxy(proxies: pd.DataFrame, proxymap: pd.DataFrame) -> pd.DataFrame:
    '''
    Of possible proxies (in proxymap) also available in outcome sumstats (extracted)
    Retain the best single proxy variant
    '''
    data = proxies.merge(proxymap, left_on='rsid', right_on='rsid_b', how='inner')
    data = data.sort_values(['rsid_a','r2'])
    out = data.groupby('rsid_a').first().reset_index()
    return out

def orient_proxy_beta(proxies: pd.DataFrame) -> pd.Series:
    '''
    The proxymap indicates the rsid_b allele correlated with the rsid_a effect allele
    Use this information (in merged input) to align the proxy beta to the original variant's effect_allele
    '''
    return np.where(proxies['correlated_effect_allele']==proxies['effect_allele'], proxies['beta'], -(proxies['beta']))

def align_betas_exposure_increasing(data: pd.DataFrame) -> pd.DataFrame:
    ''' 
    Make the exposure effect allele trait-increasing (all betas positive)
    Then align all the outcome betas to the exposure-increasing allele and calculate per-variant ratio
    '''
    data['beta_y'] = np.where(data['beta_x']>=0, data['beta_y'], -(data['beta_y']))
    data['beta_x'] = np.where(data['beta_x']>=0, data['beta_x'], -(data['beta_x']))
    return data

def merge_matches_onto_inst(signals, matches):
    signals = signals[['rsid','beta','se','pval','n','entity']]
    matches = matches[['rsid','beta','se','pval','n','entity']]
    out = signals.merge(matches, on='rsid', how='inner', validate='1:1')
    out['rsid_x'] = out['rsid']
    out['rsid_y'] = out['rsid']
    out = out.drop(columns='rsid')
    return out

def merge_proxies_onto_inst(signals, proxies):
    signals = signals[['rsid','beta','se','pval','n','entity']]
    proxies = proxies[['rsid_a','rsid','beta','se','pval','n','entity']]
    out = signals.merge(proxies, left_on='rsid', right_on='rsid_a', how='inner', validate='1:1')
    out = out.drop(columns='rsid_a')
    return out

def get_analytic_dataframe(signals: pd.DataFrame, 
                     proxymap: pd.DataFrame, 
                     ypath: str) -> pd.DataFrame:
    '''
    Wrapper function
    '''
                         
    siglist = list(zip(signals['chr'], signals['rsid']))
        
    with pd.HDFStore(ypath, 'r') as ystore:
        matches = extract_variants(siglist, ystore)
        needprox = [x for x in siglist if not x[1] in matches['rsid'].tolist()]
        possprox = find_proxies(needprox, proxymap)
        proxlist = list(zip(possprox['chr'], possprox['rsid_b']))
        proxies = keep_best_proxy(extract_variants(proxlist, ystore), proxymap)
        proxies['beta'] = orient_proxy_beta(proxies)
        
    out = pd.concat([merge_matches_onto_inst(signals, matches), merge_proxies_onto_inst(signals, proxies)])     
    return align_betas_exposure_increasing(out)
        
    
        