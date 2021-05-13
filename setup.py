from setuptools import setup

setup(
   name='Environmental_PAH_Mutagenicity',
   version='1.0',
   description='small hydrocarbon QSAR for direct acting Ames test mutagens',
   author='Trevor Sleight',
   author_email='twsleight@gmail.com',
   packages=['Environmental_PAH_Mutagenicity'],  #same as name
   install_requires=['numpy~=1.19',
                     'pandas>=1.1',
                     'xlsxwriter~=1.3',
                     'scikit-learn~=0.23',
                     'statsmodels~=0.12',
                     'openbabel>=3.1',
                     'padelpy~=0.1', 
                     
                     ], #external packages as dependencies
)
