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
                     'scikit-learn>=0.23.2',
                     'statsmodels==0.12.0'], #external packages as dependencies
)
