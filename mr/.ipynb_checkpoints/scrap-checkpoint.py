'''
Functions which are currently removed
'''

def calc_first_order_weights(se_y: pd.Series, beta_x: pd.Series) -> pd.Series:
    ''' 
    Calculate first-order (IVW) weights as used by downstream analyses (mainly radial)
    Variously denoted as (beta_x**2)/(se_y**2) or (se_y/beta_x)**-2.  These are equivalent
    '''
    # Alternate approach for first-order weights
    # Code is currently hard-coded to use inverse variance (of yvar) throughout
    return (se_y/beta_x)**-2

def calc_penalised_ivw(ratio: pd.Series,
                       ivw_weight: pd.Series,
                       betaivw: float) -> pd.Series:
    ''' Calculate the penalised weight '''
    penalty = (scipy.stats.chi2.sf(ivw_weight *(ratio - betaivw) ** 2, df=1)) * 20
    return ivw_weight * np.minimum(penalty, pd.Series(1, range(len(penalty))))