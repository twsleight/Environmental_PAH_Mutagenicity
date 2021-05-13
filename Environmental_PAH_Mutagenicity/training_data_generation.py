


#scripts to refine the data from the full CCRIS .xml file 
#to just the rows out of it that are desired. In this case just the direct
#acting mutagenic compounds. 





get_ipython().magic('reset -sf')
import pandas as pd
from Environmental_PAH_Mutagenicity.read_ccris_data import CIRconvert
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

import os

from rdkit import Chem


#This takes a few minutes
# df = convert_xml_xlsx('C:\ResearchWorkingDirectory\MutagenicityTrainingData\ccris.xml')


#if you've already written it to an excel file from a previous run use that first. 
ccrisData = pd.read_excel(r"C:\ResearchWorkingDirectory\MutagenicityTrainingData\CCRISdata.xlsx", sheet_name = 'Sheet1')

data = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\mutagenicity_data.xlsx", sheet_name = 'Sheet1')

#keep this for now just for evaluation
dataRef = data[['SMILES', 'name', 'CAS', 'result']]


casNumbers = list(data['CAS'])

#select whatever strains you want to use here. 
newAmes = ccrisData.loc[ccrisData['Strain'].str.contains('TA98|TA100')]

newAmes2 = newAmes.loc[newAmes['CAS'].isin(casNumbers)]
newAmes2.drop(labels = ['Strain', 'method', 'Unnamed: 0'], axis = 1,inplace = True)

#don't use a subset yet because therea re multiple entries
newAmes2.drop_duplicates(inplace = True)

# #need to join in the smiles
# newAmes2 = pd.merge(left=dataRef, right=newAmes2, left_on='CAS', right_on='CAS', copy = False)
# newAmes2.rename(columns={'SMILES_x':'SMILES'}, inplace = True)
# newAmes2.rename(columns={'result_x':'result'}, inplace = True)

#replace the postive and negatives with 1s and zeros

newAmes2.replace({'POSITIVE':1, 'NEGATIVE':0}, inplace = True)
newAmes2.replace({'POSITIVE; 97% PURE':1, 'NEGATIVE; TOXICITY AT 50 UG/PLATE':0}, inplace = True)
newAmes2.replace({'POSITIVE; 99.4% PURE':1, 'NEGATIVE; TOXICITY AT 100 UG/PLATE':0}, inplace = True)
newAmes2.replace({'NEGATIVE; 100% PURE':0}, inplace = True)
newAmes2.replace({'NEGATIVE; 99.8% PURE':0}, inplace = True)

newAmes2.drop_duplicates(inplace = True)


for i, row in newAmes2.iterrows():
    smiles = CIRconvert(row['CAS'])
    newAmes2.at[i,'SMILES'] = smiles
    if i < 20:
        break
    print(i)

#select just the data desired
validationStr = set(r' CcOo=()-[]123456789+@#\/H')
#molecular weight, specific strains, and number of rings. No S9 activation

finalDF2 = pd.DataFrame(columns = ['CAS', 'SMILES', 'name', 'result',  'mw'  ])


for i, row in newAmes2.iterrows():
       
    #deal with the get help's later

    curSmiles = row['SMILES']
    
    if curSmiles == 'Did not work':
        continue
    
    mol = Chem.rdmolfiles.MolFromSmiles(curSmiles)
    curSmiles = Chem.rdmolfiles.MolToSmiles(mol)
    newAmes2.at[i, 'SMILES'] = curSmiles 
    row['SMILES'] = curSmiles
   
    setSmiles = set(curSmiles)
    if set(setSmiles).issubset(validationStr):
        mol = Chem.rdmolfiles.MolFromSmiles(curSmiles)
        if mol is None:
            print(curSmiles, 'had validation errors')
            continue
        numRings = Chem.rdMolDescriptors.CalcNumRings(mol)
        aromRings = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
        mw = Chem.Descriptors.ExactMolWt(mol)    
        
        if aromRings >= 1 and numRings <= 5 and mw < 500:
    
            row['mw'] = mw
            
            #each time, check against the original dataframe to see
            #if there are other entries. 
    
            curSmiles = Chem.rdmolfiles.MolToSmiles(mol)
            tempDF = newAmes2.loc[newAmes2['SMILES'] == curSmiles]
            allResult = int(bool(sum(tempDF['result'])))
            row['result'] = allResult
            
            finalDF2 = finalDF2.append(row)
        else:
            print('dumped', curSmiles)
    else:
        print('dumped ', curSmiles)

       
             
            
#some compounds have different names or multiple CAS values, so ignore those. 
finalDF2.drop_duplicates(subset = ['SMILES', 'result'], inplace = True)
#if desired, optimize the structures in Gaussian

final_Data = finalDF2

final_Data.rename(columns={'result':'result_multi'}, inplace = True)

# writer = pd.ExcelWriter('selected_data.xlsx')
# final_Data.to_excel(writer,'Sheet1')
# writer.save()                
# writer.close()

#once the gaussian .out files have been run, use this section to extract the
#and get the padel descriptors based on the optomized structures







