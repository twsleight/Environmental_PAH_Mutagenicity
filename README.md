# Environmental_PAH_Mutagenicity   [![Build Status](https://travis-ci.org/twsleight/Environmental_PAH_Mutagenicity.svg?branch=master)](https://travis-ci.org/twsleight/Environmental_PAH_Mutagenicity)  [![Documentation Status](https://readthedocs.org/projects/environmental-pah-mutagenicity/badge/?version=latest)](https://environmental-pah-mutagenicity.readthedocs.io/en/latest/?badge=latest)

QSAR for predicting the mutagenicity of environmental metabolites of Polycyclic Aromatic Hydrocarbons, based on training data from the Ames Test TA98 and TA100 strains.

This QSAR is based on data from the Chemical Carcinogenesis Research Information System. The data can be downloaded as an xml file. 
https://www.nlm.nih.gov/databases/download/ccris.html

note: a small subset of this data is used as a test for the read xml function. Users should download the full file from this original link. 

read_ccris_data.py reads the data out of ccris, converts the cas numbers to smiles codes, and writes the data off to an excel file. CAS numbers are converted to SMILES codes by http://cactus.nci.nih.gov/chemical/structure/'

padel_Data will look in the directory it was launched from for .out files. The smiles code needs to be in the header of the file as a comment. 

if gaussian ".out" files are available, padel_Data.py will attempt to use these for the structure. 

recursive featureSelection will iterate over the specified range to select the ideal number of features, based on the specified input metrics. This is a time consuming process, so it is recommended that this be run on a cluster if possible. 

