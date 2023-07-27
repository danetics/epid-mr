from setuptools import setup

setup(
    name='epid-mr',
    version='0.1',    
    description='Tools to run python-implemented MR models using standardised HDF input',
    url='https://github.com/danetics/epid-mr',
    author='Daniel Wright',
    install_requires=['pandas', 'numpy', 'h5py','statsmodels','robustats','scipy'],
)