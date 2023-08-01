'''
Generate display labels for plotted model regression lines
'''
from .reference import modelnames

def generate_legend_label(res, model:str) -> str:
    '''
    Given a regression results objects from statsmodels
    Return a pretty string including the model's beta, SE and pvalue (plus pINT if Egger)
    This string is intended as a legend label for the plotted regression line
    '''
    # Names of fitted terms can vary, so identify the term fitted
    term = res.params.index[0]
    elements = {
        'model': f"{modelnames[model]}:",
        'beta': f"{round(res.params[term], 2)}",
        #'bse': f"({round(res.bse[term], 3)})",
        'pval': f"p={res.pvalues[term]:.1e}"
    }
    if 'Intercept' in res.params:
        elements['pval_int'] = f"pINT={res.pvalues['Intercept']:.2e}"
    return ' '.join(elements.values())

def generate_legend_label_wm(res_wm: dict) -> str:
    '''
    '''
    elements = {
        'model': f"{modelnames['wm']}:",
        'beta': f"{round(res_wm['beta'], 2)}",
        #'bse': f"({round(res_wm['se'], 3)})",
        'pval': f"p={res_wm['pval']:.1e}"
    }
    return ' '.join(elements.values())
