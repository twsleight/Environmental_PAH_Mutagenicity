
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 17:35:13 2020

@author: twsle
"""


import os
import pandas as pd
import tempfile
import shutil
from openbabel import openbabel

from padelpy import from_mdl
from padelpy import padeldescriptor
padeldescriptor(maxruntime=50000)

#*********************************************************
#this script will take the optimized structures from 
#a gaussian run and use them to calculat the topological descriptors with PADEL

        
def get_padel_data(filepath):        
        
    #copy all the files in the path directory to tmpdir
    temp_dir = tempfile.gettempdir()        
    temp_path = os.path.join(temp_dir, 'tempfile')
    
    shutil.copy2(filepath, temp_path)        
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("out", "mdl")

    #write the .mdl equivalent structure of the .out file
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, filepath)   # Open Babel will uncompress automatically       
    
    mdlOut =  'temp.mdl'    
    
    obConversion.WriteFile(mol, 'temp.mdl'  )

    obConversion2 = openbabel.OBConversion()
    obConversion2.SetInAndOutFormats("out", "can")

    obConversion2.ReadFile(mol, temp_path)   # Open Babel will uncompress automatically       

    #write the smiles off
    temp_smile_path = os.path.join(temp_dir,'temp.smi')
    obConversion2.WriteFile(mol, os.path.join(temp_dir,temp_smile_path ))
    obConversion2.CloseOutFile()
    
    curSmiles = ' '
    
    prevLine = []
    #get the smiles from the header also
    with open(filepath)as searchfile:
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
                     
            prevLine2 = prevLine             
            prevLine = line  
    
    
    #read it back in. no better way to do this. 
    a_file = open(temp_smile_path , "r")
    file_lines = a_file.readlines()
    a_file.close()
    smiles = file_lines[0].split()[0]
    # print(smiles)
    # print(curSmiles)
    
    # smiles = repr(mol.write("smi"))
    
    # #read the mdl file back in a clean it up. 
    # a_file = open(mdlOut, "r")
    # list_of_lines = a_file.readlines()
    # a_file.close()
         
    # list_of_lines[0] = "\n"
    # #put the smiles codes here
    # list_of_lines[1] = "\n"
    # list_of_lines[2] = "\n"
    
    # a_file = open(mdlOut, "w")
    # a_file.writelines(list_of_lines)
    # a_file.close()  
    # #mdl file is now cleaned up

    # #call padel and get the descriptors
    [descriptors] = from_mdl(mdlOut, timeout = 120)    
    
    #tempslice is a dataframe that contains all the descriptors from padel       
    tempSlice = pd.DataFrame.from_dict(descriptors, orient = 'index')
    # i will become the index
    
    
    # tempSlice.rename(columns={0: i}, inplace = True)
    tempSlice = tempSlice.T
    short_filename = filepath.split(r'\\')[-1]

    tempSlice['filename'] = short_filename

    # tempSlice['SMILES'] = smiles
    tempSlice['SMILES'] = curSmiles
    bigDataFrame = tempSlice

    # bigDataFrame = bigDataFrame.append(tempSlice)
    # i+=1
    # print(i)
        
    #adjust the order of the columns so that the smiles codes
    #are on the left    
    endcols = list(bigDataFrame.columns)
    # endcols.remove('newSMILES')
    endcols.remove('SMILES')

    columnList = ['SMILES'] + endcols
    bigDataFrame = bigDataFrame[ columnList]

    return bigDataFrame
