# Environmental_PAH_Mutagenicity  

QSAR for predicting the mutagenicity of environmental metabolites of Polycyclic Aromatic Hydrocarbons, based on training data from the Ames Test TA98 and TA100 strains.

This QSAR is based on data from the Chemical Carcinogenesis Research Information System. The data can be downloaded as an xml file. 
https://www.nlm.nih.gov/databases/download/ccris.html

note: a small subset of this data is used as a test for the read xml function. Users should download the full file from this original link. 

read_ccris_data.py reads the data out of ccris, converts the cas numbers to smiles codes, and writes the data off to an excel file. CAS numbers are converted to SMILES codes by http://cactus.nci.nih.gov/chemical/structure/'

padel_Data will look in the directory it was launched from for .out files. The smiles code needs to be in the header of the file as a comment. 

A test dataset is provided as the manusript for this analysis is currently in review

You will also need openbabel and rdkit installed

Note: This research is being actively used in other projects, so pushes may be made periodically but will not affect the overall functionality of the system. 

