


#main file to run the scripts that generate the data based off of the 

#guassian .out files

#and Sauron forged in secret a master ring...

from IPython import get_ipython;   
get_ipython().magic('reset -sf')

import pandas as pd
import glob
import os

from rdkit import Chem

from Environmental_PAH_Mutagenicity.read_ccris_data import convert_xml_xlsx
# from Environmental_PAH_Mutagenicity.read_ccris_data import CIRconvert
from Environmental_PAH_Mutagenicity.gethlenergy import gethlenergy
from Environmental_PAH_Mutagenicity.get_padel_data import get_padel_data
# from Environmental_PAH_Mutagenicity.gethlenergy import gethlenergy
from Environmental_PAH_Mutagenicity.read_IP_EA_functions import read_neut_energy
from Environmental_PAH_Mutagenicity.read_IP_EA_functions import read_posneg_energy


#a handy script for creating the guassian inputs from smiles codes
#not used in this file
def initalize_guass_input(filename, smiles):
    
    m = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m)
    
    AllChem.EmbedMolecule(m,randomSeed=0xf00d)  
    
    # AllChem.MMFFOptimizeMolecule(m)
    
    coord_3d = Chem.MolToMolBlock(m)

    coordArray = coord_3d.splitlines()
    

    actualCoords = []
    for c in coordArray:
        if len(c) == 69:
            chunks = c.split()
            actualCoords.extend([chunks[3],'          ',chunks[0], '          ',chunks[1], '          ',chunks[2], '\n'])
            

    #I want to overwrite the entire file, so I need to open it back up in write mode. 
    list_of_lines = list(range(0,8))
    list_of_lines[0] = "%Nproc=6\n"
    list_of_lines[1] = "%mem=20GB\n"
    list_of_lines[2] = "#b3lyp/6-311g(d,p) opt=maxcycle=640 freq nosymm scf=(maxcycle=999) gfinput iop(6/7=3) "

    list_of_lines[3] = '\n'     
    list_of_lines[4] = '\n'
    list_of_lines.insert(5, (smiles +  "\n"))
    
    list_of_lines[6] = "\n0 1\n"  
    list_of_lines[7] = "".join(actualCoords)
    list_of_lines[8] = '\n \n \n'

    a_file = open(filename, "w")
    a_file.writelines(list_of_lines)
    a_file.close()

       
    return


# testDF = convert_xml_xlsx("test_CCRIS_data.xml")

runningList = []
bigDataFrame = pd.DataFrame()

path = r'C:\ResearchWorkingDirectory\MutagenicityDataMar\completed'
verticals_path = r'C:\ResearchWorkingDirectory\MutagenicityDataMar\posneg'


i = 0
for filename in glob.glob(os.path.join(path,'*.OUT')): 
    #get the smiles from the file, if we already have this smiles, skip and go to next.
    
    #Use this to get the appropriate verticals
    short_filename = filename.split('\\')[-1]
    base_filename = short_filename.split('.')[0]    
    
    if os.path.exists(os.path.join(verticals_path, (base_filename +'_pos.out') )):
    
        ip_eng = read_posneg_energy(verticals_path, (base_filename +'_pos.out'))
        
    else:
         print((base_filename +'_pos.out') + ' not found')
    
    if os.path.exists(os.path.join(verticals_path, (base_filename +'_neg.out') )):        
    
        ea_eng = read_posneg_energy(verticals_path, (base_filename +'_neg.out'))
    else:
        print((base_filename +'_neg.out') + ' not found')


    prevLine = []
    curSmiles = []    
    with open(filename)as searchfile:
        i = i+1
        print(i)
        #************************  
        p = 'not found'
        c = 'not found'
        t = 'not found'
        for line in searchfile:
        #***************************************************************
            #Get the smiles code that we have in the top of every input card
            left,sep,right=line.partition('Symbolic Z-matrix:') 
            
            #as long as the smiles code is in the line this code will find it. 
            #there can be other stuff there also, as long as the smiles code
            #is seperated by at least one space. 
            if sep:
    #                print(sep)
    #                    print(prevLine)
                p = prevLine2. split()[-1]
                
                if set(p).issubset(r' CcOo=()-[]123456789+@#H\//') and ('c' in set(p) or 'C' in set(p) ):  
#                            print(p)
                    curSmiles = p
                    # print(prevLine)
                else:
                    print('issue with ', p, 'filename ', filename)
                     
            prevLine2 = prevLine             
            prevLine = line       
   
    #check curSmiles 
    mol = Chem.rdmolfiles.MolFromSmiles(curSmiles)
    curSmiles = Chem.rdmolfiles.MolToSmiles(mol)    
    
    
    # if curSmiles in runningList:
    #     continue
    
    paData = get_padel_data(filename)
    
    #Use this to get the appropriate verticals
    short_filename = filename.split('\\')[-1]
    base_filename = short_filename.split('.')[0]
    
    ip_eng = 0
    ea_eng = 0

    
    if os.path.exists(os.path.join(verticals_path, (base_filename +'_pos.out') )):
    
        ip_eng = read_posneg_energy(verticals_path, (base_filename +'_pos.out'))
        
    else:
         print((base_filename +'_pos.out') + ' not found')
    
    if os.path.exists(os.path.join(verticals_path, (base_filename +'_neg.out') )):        
    
        ea_eng = read_posneg_energy(verticals_path, (base_filename +'_neg.out'))
    else:
        print((base_filename +'_neg.out') + ' not found')



    HLdata = gethlenergy(path, short_filename)
    
    paData['HOMO'] = HLdata[0]
    paData['LUMO'] = HLdata[1]

    #add the IP and EA data right here
    neutral_eng = read_neut_energy(path, filename)
    paData['Neutral'] = neutral_eng
    
    paData['IP'] = ip_eng


    paData['EA'] = ea_eng

    
    bigDataFrame = bigDataFrame.append(paData)
    
        
    runningList.append(curSmiles)
    #merge this into the existnig data frame based on the smiles code


    # short_filename = filename.split(r'\\')[-1]
  
    # idx = final_Data.index[final_Data['SMILES']==curSmiles].tolist()
    # mw = Chem.Descriptors.ExactMolWt(mol)
    # if len(idx)>0:
    #     final_Data.at[idx[0], 'HOMO'] = HLdata[0]
    #     final_Data.at[idx[0], 'LUMO'] = HLdata[1]
    #     final_Data.at[idx[0], 'filename'] = short_filename
        
        
    # else:
    #     pass
    #     print(curSmiles)
    #     print(filename.split(r'\\')[-1])
    #     # print('molecular weight is ', mw)

bigDataFrame['HLgap'] = (bigDataFrame['LUMO']-bigDataFrame['HOMO'])*27.2114

bigDataFrame['IP_eV'] = (bigDataFrame['Neutral'] - bigDataFrame['IP'])*27.2114

bigDataFrame['EA_eV'] = (bigDataFrame['Neutral'] - bigDataFrame['EA'])*27.2114


endcols = list(bigDataFrame.columns)
endcols.remove('SMILES')
# endcols.remove('result')
endcols.remove('HOMO')
endcols.remove('LUMO')
endcols.remove('IP_eV')
endcols.remove('EA_eV')
endcols.remove('HLgap')
# endcols.remove('result')


#move some columsn to the front
columnList = ['SMILES','HOMO', 'LUMO', 'IP_eV', 'EA_eV', 'HLgap'] + endcols
bigDataFrame2 = bigDataFrame[ columnList]


#this the reference data. Use this to get name and CAS. CAS is very useful later
# bonus_data = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\beforeMarchReset\TA98TA100_RawData3.xlsx", sheet_name = 'Sheet1')


bonus_data = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\mutagenicity_data.xlsx", sheet_name = 'Sheet1')



final_results = bonus_data[['SMILES', 'name', 'CAS', 'result']]

# bigDataFrame2 = current_data.merge(final_results, on='SMILES', how='inner')



bigDataFrame2 = bigDataFrame2.merge(final_results, on='SMILES', how='inner')

endcols = list(bigDataFrame2.columns)
endcols.remove('Neutral')
endcols.remove('result')
endcols.remove('IP')
endcols.remove('EA')
endcols.remove('filename')
endcols.remove('CAS')
endcols.remove('name')

columnList = ['Neutral','IP', 'EA', 'filename', 'CAS', 'name', 'result'] + endcols
bigDataFrame2 = bigDataFrame2[ columnList]




#Write off the final data

# writer = pd.ExcelWriter('b3lyp_data_2.xlsx')
# bigDataFrame2.to_excel(writer,'Sheet1')
# writer.save()                
# writer.close()

