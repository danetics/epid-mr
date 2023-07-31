'''
Generate text to be included in plot area
'''

from statsmodels.stats import weightstats

def get_variant_counts(data: pd.DataFrame, highlight: pd.Series) -> str:
    '''
    Given an analytic dataframe, return the counts of variants therein as a formatted string for display on the plot
    If a highlight pd.Series is also passed, remove those variants from the data before counting
    '''
    # Use case for removing the highlighted variants is the filtered models
    data = data[~highlight]
    return f"Variants={len(data)} ({len(data[data['rsid_x']!=data['rsid_y']])} proxies)"

def get_qstats(data, reglines: list):
    if reglines==['ivw_radial_filtered','egger_radial_filtered']:
        cochran = 
        rucker = 
        return f"{cochran}\n{rucker}"
    res = weightstats.wls_meta(data['ratio'], data['ratio_se']**-2)
    return f"Q={} p={} (I)"
        

def add_textbox(data: pd.DataFrame, res, reglines: list):
    '''
    Generate a textbox including variant counts and heterogeneity statistics for a given model
    '''
    counts = f"Variants = {res[model]['reg'].nobs} ({len(data[data['rsid_x']!=data['rsid_y']])} proxies)"
    if model in ['radial','egger_radial']:
        cochran = f"Cochran's Q = {data['rad_cochranq'].sum()} p={chi2.sf(data['rad_cochranq'].sum(), len(data)-1):.2e}"
        rucker = f"Rucker's Q = {data['rad_ruckerq'].sum()} p={chi2.sf(data['rad_ruckerq'].sum(), len(data)-1):.2e}"
        return f"{counts}\n{cochran}\n{rucker}"
    cochran = 