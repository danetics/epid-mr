'''
Self-explanatory subfunctions for handling HDF files
'''

import h5py
import pandas as pd

def read_metadata(path: str) -> dict:
    with h5py.File(path, 'r') as file:
        return dict(file.attrs.items())

def load_signals(path: str) -> pd.DataFrame:
    with pd.HDFStore(path, 'r') as x:
        return x.get('signals')

def extract_instrument(path: str, 
                       instrument: pd.DataFrame) -> pd.DataFrame:
    ''' Extract a given instrument from a set of HDF sumstats, returning a merged analytic dataframe '''
    data = pd.DataFrame()
    with pd.HDFStore(path, 'r') as store:
        for k,v in instrument.groupby('chr')['bp'].apply(list).to_dict().items():
            lookup = store.select(f'chr{k}', where="bp = v")
            data = pd.concat([data, lookup])
        
    return instrument.merge(data, on=['bp'], how='inner')
    # Default suffixes are [x,y], which fits the x and y vars based on above merge syntax
    # Will produce beta_x, se_x etc, which is desired naming