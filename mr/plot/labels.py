'''
Generate display labels for plotted model regression lines
'''

# Display names for each model


def generate_legend_label(res: sm.regression.linear_model.RegressionResultsWrapper, model:str) -> str:
    '''
    Given a regression results objects from statsmodels
    Return a pretty string including the model's beta, SE and pvalue (plus pINT if Egger)
    This string is intended as a legend label for the plotted regression line
    '''
    # Names of fitted terms can vary, so identify the term fitted
    term = res.params.index[0]
    elements = {
        'model' = f"{modelnames[model]}:",
        'beta' = f"{round(res.params[term], 2)}",
        'bse' = f"({round(res.bse[term], 3)})",
        'pval' = f"p={res.pvalues[term]:.2e}"
    }
    if 'Intercept' in res.params:
        elements['pval_int'] = f"pINT={res.pvalues['Intercept']:.2e}"
    return ' '.join(elements.values())
