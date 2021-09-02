"""
Created on Thu Sep 17 18:40:47 2020

@author: twsle
"""


# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.lines import Line2D
# import matplotlib.pyplot as plt
import matplotlib.colors

def getNonZeroColumns(data):
    headers = list(data.columns)
    dropnondata = []
    for h in headers: 
        if data[h].dtype != np.number:
             dropnondata.append(h)
        elif sum(data[h]) == 0:
            dropnondata.append(h)           
    headers = list(set(headers) - set(dropnondata))    
        
    return headers

# Load in the data
data = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\mutagenicity_data.xlsx", sheet_name = 'Sheet1')
gramData = pd.read_csv(r"C:\ResearchWorkingDirectory\extrapapersData\grammatica\grammaticaDatabig.csv")
pgData = pd.read_csv(r"C:\ResearchWorkingDirectory\extrapapersData\perezgarrido\perezgarridoDatabig.csv")
haoData = pd.read_csv(r"C:\ResearchWorkingDirectory\extrapapersData\hao\haoDatabig2.csv")
reenuData = pd.read_csv(r"C:\ResearchWorkingDirectory\extrapapersData\Reenu\reenuDatabig2.csv")
debnathData =  pd.read_csv(r"C:\ResearchWorkingDirectory\extrapapersData\debnath\debnath_Databig2.csv")
wangData =  pd.read_csv(r"C:\ResearchWorkingDirectory\extrapapersData\wang\wang_Databig2.csv")
kimData =  pd.read_csv(r"C:\ResearchWorkingDirectory\extrapapersData\kim\kim_Databig2.csv")

headers = set(data.columns)
headershao = set(haoData.columns)
headerspg = set(pgData.columns)
headersG = set(gramData.columns)
headersreenu = set(reenuData.columns)
headersdebnath = set(debnathData.columns)
headerswang = set(wangData.columns)
headerskim = set(kimData.columns)

finalheaders = list(headers.intersection(headersG).
                    intersection(headerspg).
                    intersection(headershao).
                    intersection(headersreenu).
                    intersection(headersdebnath).
                    intersection(headerswang).
                    intersection(headerskim))   

data0 = data[finalheaders]
data0['source'] = ['#CC3311']* len(data0) #data from this study


#data from other studies
gramColor = '#808080'
gramData = gramData[finalheaders]
gramData['source'] = [gramColor]* len(gramData)


pgColor = '#808080'
pgData = pgData[finalheaders]
pgData['source'] = [pgColor]* len(pgData)

haoColor = '#808080'
haoData = haoData[finalheaders]
haoData['source'] = [haoColor]* len(haoData)

reenuColor = '#808080'
reenuData = reenuData[finalheaders]
reenuData['source'] = [reenuColor]* len(reenuData)

debnathColor = '#808080'
debnathData= debnathData[finalheaders]
debnathData['source'] = [debnathColor]* len(debnathData)

wangColor = '#808080'
wangData= wangData[finalheaders]
wangData['source'] = [wangColor]* len(wangData)

kimColor = '#808080'
kimData= kimData[finalheaders]
kimData['source'] = [kimColor]* len(kimData)


#put the two bigger studies in front, then the smaller studies, my data at the back. 
finaldata = pd.concat([pgData, debnathData, wangData, kimData, haoData,reenuData, gramData,data0])
finaldata = finaldata.dropna(axis = 1)

#Get the colors for the clusters laster
finalheaders.append('source')
colorsSource = finaldata['source']

finaldata = finaldata.drop(columns = ['source'])

headers2 = getNonZeroColumns(finaldata)
X = StandardScaler().fit_transform(finaldata[headers2])

# Create a PCA instance: pca\
pca = PCA(n_components=10)
principalComponents = pca.fit_transform(X)
# Plot the explained variances
features = range(pca.n_components_)
plt.bar(features, pca.explained_variance_ratio_, color='black')
plt.xlabel('PCA features')
plt.ylabel('variance %')
plt.xticks(features)

# Save components to a DataFrame, these are the
# principle components of everything
PCA_components = pd.DataFrame(principalComponents)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(PCA_components[0], PCA_components[1], PCA_components[2], alpha=.33, c = colorsSource)
plt.xlabel('PCA 1')
plt.ylabel('PCA 2')


allalpha = 0

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

#Perez-Garrido
rgbCol = matplotlib.colors.to_rgb(pgColor)
ax.scatter(PCA_components[1][0:219], PCA_components[2][0:219], 
           color = (rgbCol[0],rgbCol[1],rgbCol[2],allalpha), 
           edgecolor = (0,0,0,1),linewidth = 0.75, marker = 'v', zorder = 0)


#Debnath
rgbCol = matplotlib.colors.to_rgb(debnathColor)
ax.scatter(PCA_components[1][220:380], PCA_components[2][220:380],  
           color = (rgbCol[0],rgbCol[1],rgbCol[2],allalpha), 
           edgecolor = (0,0,0,1),linewidth = 0.75,marker = '^', zorder = 0)

#Wang
rgbCol = matplotlib.colors.to_rgb(wangColor)
ax.scatter(PCA_components[1][381:407], PCA_components[2][381:407],  
          color = (rgbCol[0],rgbCol[1],rgbCol[2],allalpha),
           edgecolor = (0,0,0,1), linewidth = 0.75,marker = '>', zorder = 2)

#Kim
rgbCol = matplotlib.colors.to_rgb(kimColor)
ax.scatter(PCA_components[1][408:438], PCA_components[2][408:438],  
          color = (rgbCol[0],rgbCol[1],rgbCol[2],allalpha), 
          edgecolor = (0,0,0,1),linewidth = 0.75, marker = '<', zorder = 3)

#Hao
rgbCol = matplotlib.colors.to_rgb(haoColor)
ax.scatter(PCA_components[1][439:487], PCA_components[2][439:487],  
            color = (rgbCol[0],rgbCol[1],rgbCol[2],allalpha), 
            edgecolor = (0,0,0,1),linewidth = 0.75, marker = 'P', zorder = 3)

#Reenu
rgbCol = matplotlib.colors.to_rgb(reenuColor)
ax.scatter(PCA_components[1][488:539], PCA_components[2][488:539],  
            color = (rgbCol[0],rgbCol[1],rgbCol[2],allalpha), 
            edgecolor = (0,0,0,1),linewidth = 0.75, marker = 's', zorder =2)

#Grammatica
rgbCol = matplotlib.colors.to_rgb(gramColor)
ax.scatter(PCA_components[1][540:565], PCA_components[2][540:565],  
            color = (rgbCol[0],rgbCol[1],rgbCol[2],allalpha),
            edgecolor = (0,0,0,1),linewidth = 0.75, marker = 'D', zorder = 2)


#plot the current study data on top
rgbCol = matplotlib.colors.to_rgb('#CC3311')
ax.scatter(PCA_components[1][565:], PCA_components[2][565:],  
            color = (rgbCol[0],rgbCol[1],rgbCol[2],.9),
            edgecolor = (0,0,0,1),linewidth = 0.75, marker = 'o', zorder = 1)



plt.xlabel('Principal Component 1', fontsize = 16)
plt.ylabel('Principal Component 2', fontsize = 16)

custom_lines = [Line2D([0], [0], color='#CC3311',linewidth = 1, markersize = 10, markeredgecolor = 'black', marker = 'o',lw = 0, alpha = 0.9),  #this study
                Line2D([0], [0], color='#FFFFFF',linewidth = 1, markersize = 10, markeredgecolor = 'black', marker = 'P',lw=0, alpha = 1),  #Hao
                Line2D([0], [0], color='#FFFFFF',linewidth = 1, markersize = 10, markeredgecolor = 'black', marker = 's',lw=0, alpha = 1), #Reenu
                Line2D([0], [0], color='#FFFFFF',linewidth = 1,  markersize = 10, markeredgecolor = 'black', marker = 'v',lw=0, alpha = 1), #Perez-Garrido
                Line2D([0], [0], color='#FFFFFF',linewidth = 1, markersize = 10, markeredgecolor = 'black', marker = '>', lw=0, alpha = 1), #Wang
                 Line2D([0], [0], color='#FFFFFF',linewidth = 1, markersize = 10, markeredgecolor = 'black',  marker = 'D',lw=0, alpha = 1), #Grammatica              
                 Line2D([0], [0], color='#FFFFFF',linewidth = 1, markersize = 10, markeredgecolor = 'black', marker = '<',lw=0, alpha = 1), #Kim        
                Line2D([0], [0], color='#FFFFFF',linewidth = 1, markersize = 10, markeredgecolor = 'black',marker = '^',lw=0, alpha = 1)]  #Debnath


#put these in chronological order
ax.legend(custom_lines, ( 'This Study', 
                         'Hao, 2019', 
                         'Reenu, 2014',
                         'Perez-Garrido, 2010',                       
                         'Wang, 2009',                         
                         'Grammatica, 2007 ',
                         'Kim, 2006',                       
                         'Debnath, 1992' ))

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.set_xlim([-21,50])
ax.set_ylim([-21,37])
ax.spines['top'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax = plt.axes()
ax.grid(False)
plt.show()

fig.savefig('otherQSARS_bw.png', format='png', dpi=1200)

