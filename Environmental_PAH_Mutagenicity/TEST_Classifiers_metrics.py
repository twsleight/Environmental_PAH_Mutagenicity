

import pandas as pd
from sklearn.metrics import confusion_matrix 
from rdkit import Chem


def standardize_SMILES(df, column_smiles):
     
    for i, row in df.iterrows():
        curSmiles = row[column_smiles]
        mol = Chem.rdmolfiles.MolFromSmiles(curSmiles)
        curSmiles2 = Chem.rdmolfiles.MolToSmiles(mol)
        
        df.at[i, column_smiles] = curSmiles2
    return df



df = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\TEST_predictions.xlsx", sheet_name = 'Sheet1')
df['Pred_Result'] = (df['Pred_Consensus'] >= 0.5).astype(int)



conf = confusion_matrix(df['result'], df['Pred_Result'])
tn = conf[0,0]
tp = conf[1,1]
fn = conf[1,0]
fp = conf[0,1]

acc = (tp+tn) / (tp + tn + fn + fp)

print('Accuracy is ', acc) 
print('Precision is ', tp/(tp+fp))
prec = tp/(tp+fp)
print('Recall/Sensitivity ', tp/(tp+fn))  
rec =  tp/(tp+fn)
print('Specificity is ', tn/(tn+fp))
print('F1 is ', 2*(prec*rec)/(prec+rec))

