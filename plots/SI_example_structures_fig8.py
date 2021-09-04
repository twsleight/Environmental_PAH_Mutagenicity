# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 09:32:14 2020

@author: twsle
"""

#read in the Reference data and sort it into positives and negatives
#**************************************************************

import os
import pandas as pd
import matplotlib.pyplot as plt
import time
from rdkit import Chem
from rdkit.Chem.Draw import DrawingOptions


#Draing options for the display of the MCS structures at the top of the file
DrawingOptions.atomLabelFontSize = 120
DrawingOptions.dotsPerAngstrom = 120
DrawingOptions.bondLineWidth = 16.0     #default was 20
DrawingOptions.colorBonds   =  False
DrawingOptions.noCarbonSymbols  = True
DrawingOptions.defaultColor = (1,0,0)
DrawingOptions.dblBondOffset = 0.35

# DrawingOptions.elemDict= {0: (0.5, 0.5, 0.5), 1: (0.55, 0.55, 0.55), 
#                           7: (0, 0, 1), 8: (1,0.5,0), 9: (0.2, 0.8, 0.8), 
#                           15: (1, 0.5, 0), 16: (0.8, 0.8, 0), 17: (0, 0.8, 0), 
#                           35: (0.5, 0.3, 0.1), 53: (0.63, 0.12, 0.94)}

#highlight color
hcolor = (1,0.5,0)
#set the directory
os.chdir(r'C:\ResearchWorkingDirectory\MutagenicityTrainingData')

#load the data
parent = os.path.join(os.path.abspath(__file__), os.pardir)

filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'mutagenicity_data.xlsx'))
dataOrig = pd.read_excel(filename, sheet_name = 'Sheet1')


filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'cluster_data.xlsx'))
clust_DF = pd.read_excel(filename, sheet_name = 'cluster_assignments')


orig_Data = dataOrig.loc[clust_DF['cluster']==1]
# orig_Data = dataOrig
allData = orig_Data.reset_index(drop = True)


# allData = allData.loc[allData['result']==1].reset_index(drop = True)

#small cluster non mut
#small cluster mut
#large clsuter non mut
#large cluster mut

# selectedSmiles = ['','CC(C)Cc1ccc(C(C)C(=O)O)cc1', 'COc1cc(OC)c(OC)cc1C=O','CCCCCCc1ccc(O)cc1O',
                  
#                   'Oc1cc(O)c2ccccc2c1', 'O=CC(=O)c1ccccc1', 'c1ccc2c(c1)COC2',
                  
#                   'O=C1c2ccccc2C(=O)c2c(O)cc(O)cc21','O=C1c2ccccc2-c2ccccc21','O=C(O)c1occ(-c2ccccc2)c1-c1ccccc1',
                  
#                   'CC(=O)OCc1c2ccccc2cc2ccccc12', 'COc1cc(O)c2c(O)c3c(=O)oc(C)cc3c(OC)c2c1', 'COC(=O)c1cc2c(c3c1ccc1c(OC)cccc13)OCO2',
#                   ]

colormatrix = [0, 0,0,0,1,1,1,0,0,0,1,1,1]

nrow = 3; ncol = 5;

fig, axs = plt.subplots(nrows=nrow, ncols=ncol,subplot_kw={'xticks': [], 'yticks': []})

i = 1

startpoint = 1

os.chdir(r'C:\ResearchWorkingDirectory\tempImages')

# import matplotlib.cbook as cbook
#*****************************
for ax in axs.reshape(-1):
       
    smiles = allData.loc[i+startpoint, 'SMILES']
    # smiles = selectedSmiles[i]
    
    # homor = allData.loc[allData['SMILES'] == smiles]['HOMO'].reset_index()  
    # homo = homor['HOMO'][0]*27.2114

    # spmin5r = allData.loc[allData['SMILES'] == smiles]['SpMin5_Bhi'].reset_index()  
    # spmin5 = spmin5r['SpMin5_Bhi'][0]

    hlgap = allData.loc[i+startpoint, 'HLgap']    
    rdf55 = allData.loc[i+startpoint, 'RDF55s']
    
    print(i)
    # print(homo)
    
    # mutr = allData.loc[allData['SMILES'] == smiles]['result'].reset_index()  
    # mut = mutr['result'][0]

    mut =  allData.loc[i+startpoint, 'result']
    
    # mut = colormatrix[i]
    
    if mut == 1:
        DrawingOptions.elemDict= {0: (0.5, 0.5, 0.5), 
                                  1: (0.55, 0.55, 0.55), 
                                  7: (0, 0, 1), 8: (1,0,0), 
                                  9: (0.2, 0.8, 0.8), 
                                  15: (1, 0.5, 0), 
                                  16: (0.8, 0.8, 0), 
                                  17: (0, 0.8, 0), 
                                  35: (0.5, 0.3, 0.1), 
                                  53: (0.63, 0.12, 0.94)}
        
        DrawingOptions.defaultColor = (1,0,0)
        mutColor = 'black'

    else:
        DrawingOptions.elemDict= {0: (0.5, 0.5, 0.5), 
                                  1: (0.55, 0.55, 0.55), 
                                  7: (0, 0, 1), 8: (0,0,0), 
                                  9: (0.2, 0.8, 0.8), 
                                  15: (0, 0, 0), 
                                  16: (0.8, 0.8, 0), 
                                  17: (0, 0.8, 0), 
                                  35: (0.5, 0.3, 0.1), 
                                  53: (0.63, 0.12, 0.94)}

        DrawingOptions.defaultColor = (0,0,0)
        mutColor = 'black'
   
    ax = plt.subplot(nrow,ncol,i)
     
    m = Chem.MolFromSmiles(smiles)
    
    print(smiles)
    
    filename = 'temp.png'
    
    print(filename)
    
    img = Chem.Draw.MolToFile(m,filename, size=(1200, 1200),fitImage=False)
    # with cbook.get_sample_data('temp.png') as image_file:
    image = plt.imread(filename)
    
    #pause to make sure the image is finished
    time.sleep(1)

    # homo = 'HOMO: {}'.format(round(homo, 2))+ ' eV'
    # spmin5 = 'SpMin5_Bhi: {}'.format(round(spmin5, 2)) 
    
    # plt.text(1,1400,homo, fontsize = 6, weight = 'bold')      
    # plt.text(1,1600, spmin5, fontsize = 6, weight = 'bold')      
    
    plt.imshow(image)
    
    
    
    
    ax.axis("on")
    plt.setp(ax.spines.values(), linewidth=1.25, color = mutColor)
    

    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])    
    plt.show()

    i += 1
#*****************************
plt.tight_layout()
plt.show()


# fig.savefig('SI examples other.jpg', format='jpg', dpi=1200)



