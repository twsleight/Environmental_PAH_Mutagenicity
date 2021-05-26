
import pandas as pd
from sklearn.metrics import confusion_matrix 
from rdkit import Chem
import os


def standardize_SMILES(df, column_smiles):
     
    for i, row in df.iterrows():
        curSmiles = row[column_smiles]
        mol = Chem.rdmolfiles.MolFromSmiles(curSmiles)
        curSmiles2 = Chem.rdmolfiles.MolToSmiles(mol)
        
        df.at[i, column_smiles] = curSmiles2
    return df



# #translate the mutagenicity data into numbers
# for i, row in df.iterrows():
#     mut = row['Mutagenicity (Ames test) CONSENSUS model - assessment']
    
#     df.at[i, 'Pred_Result'] = 0  
    
#     if mut.find("Mutagenic") == 0: 
#         df.at[i, 'Pred_Result'] = 1      
    
#     if mut.find("NON-Mutagenic") > 0:
#         df.at[i, 'Pred_Result'] = 0

  
parent = os.path.join(os.path.abspath(__file__), os.pardir)

filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'vega_predictions.xlsx'))

df = pd.read_excel(filename, sheet_name = 'Sheet1')
      
    
#calculate the confusion matrix
conf = confusion_matrix(df['result_from_database'], df['Pred_Result'])
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




