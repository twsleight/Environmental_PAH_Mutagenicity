import os
import pandas as pd
import numpy as np

from Environmental_PAH_Mutagenicity.classify_mutagens_functions import regress_kfold

def test_regress_kfold():
    parent = os.path.join(os.path.abspath(__file__), os.pardir)

    path = os.path.abspath(os.path.join(parent,'..', 'data', 'testdata.xlsx'))
    
    finaldata = pd.read_excel(path, sheet_name = 'descriptors')
    y = pd.read_excel(path, sheet_name = 'y')
        
    confs_results = regress_kfold(finaldata, y['result'], n_splits=3)
    
    confs = confs_results[0]
    assert (confs == np.array([[319,  67],[ 18, 113]])).all()
