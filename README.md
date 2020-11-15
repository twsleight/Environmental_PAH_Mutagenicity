# Environmental_PAH_Mutagenicity   [![Build Status](https://travis-ci.org/twsleight/Environmental_PAH_Mutagenicity.svg?branch=master)](https://travis-ci.org/twsleight/Environmental_PAH_Mutagenicity)

QSAR for predicting the mutagenicity of environmental metabolites of Polycyclic Aromatic Hydrocarbons, based on training data from the Ames Test TA98 and TA100 strains.

This QSAR is based on data from the Chemical Carcinogenesis Research Information System. The data can be downloaded as an xml file. 
https://www.nlm.nih.gov/databases/download/ccris.html

read_ccris_data.py reads the data out of ccris, converts the cas numbers to smiles codes, and writes the data off to an excel file. 

padel_Data

if gaussian ".out" files are available, padel_Data.py will attempt to use these for the structure. 



