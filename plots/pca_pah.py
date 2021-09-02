


#code to generat Figure 3 in the main text and SI figures  S1 - S6

# Imports
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


from sklearn import decomposition
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import davies_bouldin_score
from sklearn.metrics import silhouette_score


parent = os.path.join(os.path.abspath(__file__), os.pardir)

filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'mutagenicity_data.xlsx'))

data = pd.read_excel(filename, sheet_name = 'Sheet1')

data.dropna(axis = 0, inplace = True)

headers2 = list(data.columns)
headers2.remove('filename')
headers2.remove('IP')
headers2.remove('EA')
headers2.remove('Neutral')
headers2.remove('SMILES')
headers2.remove('result')
headers2.remove('CAS')
headers2.remove('name')
headers2.remove('Data Source')
data = data[headers2]

plt.rcParams['axes.grid'] = False

#remove the zero data with calcuation failures
dropnondata = []
for h in headers2: 
    if sum(data[h]) == 0:
        dropnondata.append(h)
headers2 = list(set(headers2) - set(dropnondata))

data = data[headers2]
data.reset_index(drop = True, inplace = True)

X = data[headers2]
X = StandardScaler().fit_transform(X)


#SI Figures S1
#Explained by each principle component
plt.figure(1)

# Create a PCA instance: pca
pca = PCA(n_components=10)
principalComponents = pca.fit_transform(X)

pca_decomp = decomposition.PCA(n_components=10)
pca_decomp.fit_transform(X)

# Plot the explained variances
features = range(1, (pca.n_components_+1),1)
plt.bar(features, pca.explained_variance_ratio_, color='black')
plt.xlabel('PCA features')
plt.ylabel('Variance %')
plt.xticks(features)
plt.grid(None)

# Save the principle components to a DataFrame
PCA_components = pd.DataFrame(principalComponents)
fig =  plt.gcf()
plt.grid(b=None)
# fig.savefig('S1_variance_per_PC.png', format='png', dpi=1200)

#Try up to 10 clusters
ks = range(1,10)
inertias = []
silouts = []
davies_b = []


for k in ks:
    # Create a KMeans instance with k clusters: model
    model = KMeans(n_clusters=k)

    # use the first 5 PCA components
    X_pc = PCA_components.iloc[:,:5]
    
    #apply kmeans cluster
    model.fit(X_pc)
    cluster_labels = model.fit_predict(X_pc)
    
    #ag or spectral  
    if k >1:
                        
        cur_silhout = silhouette_score(X_pc, cluster_labels)
        silouts.append(cur_silhout)
        
        cur_davies_b = davies_bouldin_score(X_pc, cluster_labels)
        davies_b.append(cur_davies_b)

    # Append the inertia to the list of inertias
    print(k)
    inertias.append(model.inertia_)
 
#SI Figure 2
plt.figure()   
plt.plot(ks, inertias, '-o', color='black')
plt.xlabel('Number of Clusters')
plt.ylabel('Within-Cluster Sum of Squres(Inertia)')
plt.xticks(ks)
plt.tight_layout()
ax = plt.axes()
ax.grid(False)
plt.show()


fig =  plt.gcf()
# fig.savefig('S2_elbow.png', format='png', dpi=1200)

#SI Figure 3
plt.figure()
plt.plot(range(2,2+len(silouts)), silouts, '-o', color='black')
plt.xlabel('Number of Clusters')
plt.ylabel('Silhouettes')
plt.xticks(ks)
ax = plt.axes()
ax.grid(False)

plt.show()
fig =  plt.gcf()
# fig.savefig('S3_silhouette.png', format='png', dpi=1200)

#SI Figure 4
plt.figure()
plt.plot(range(2,2+len(davies_b)), davies_b, '-o', color='black')
plt.xlabel('Number of Clusters')
plt.ylabel('Davies-Bouldin')
plt.xticks(ks)
ax = plt.axes()
ax.grid(False)

plt.show()
fig =  plt.gcf()
# fig.savefig('S4_daviesbouldin.png', format='png', dpi=1200)




#Loading plots for the PC analysis

loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
loading_matrix = pd.DataFrame((loadings[:,0:2]), 
                              columns=['PC1', 'PC2'], 
                              index=range(1,len(data.columns)+1))

loading_matrix['name'] = list(data.columns)


#organize the data
plot_matrix = loading_matrix[['name', 'PC1', 'PC2']]

plot_matrix.sort_values(by = ['PC1'], inplace = True, ascending = False)
plot_final = plot_matrix[['name', 'PC1', 'PC2']][0:10]

plot_matrix.sort_values(by = ['PC1'], inplace = True, ascending = True)
plot_final = pd.concat([plot_final, plot_matrix[['name', 'PC1', 'PC2']][0:10] ])

plot_matrix.sort_values(by = ['PC2'], inplace = True, ascending = True)
plot_final = pd.concat([plot_final, plot_matrix[['name', 'PC1', 'PC2']][0:10] ])

plot_matrix.sort_values(by = ['PC2'], inplace = True, ascending = False)
plot_final = pd.concat([plot_final, plot_matrix[['name', 'PC1', 'PC2']][0:10] ])


#typically no duplicates, but drop them if there are any
plot_final.drop_duplicates(inplace = True)



plt.figure( figsize = (10,10) )
axb = fig.add_subplot(111)

#create a bar chart
y_pos = np.arange(len(plot_final))

# Create horizontal bars
plt.barh(y_pos+ 0.2, plot_final['PC1'], height = 0.4, color = 'blue')
plt.barh(y_pos-0.2,  plot_final['PC2'], height = 0.4, color = 'orange')

import matplotlib.patches as mpatches

#label this
red_patch = mpatches.Patch(color='blue', label='PC1')
yellow_patch = mpatches.Patch(color='orange', label='PC2')

plt.legend(handles=[red_patch, yellow_patch])
plt.xlabel('PC Loadings', fontsize =18)

# Create names on the x-axis
plt.yticks(y_pos, plot_final['name'])
 
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=12) 

fig = plt.gcf()
# Show graphic
plt.show()
# fig.savefig('PC_loadings.png', format='png', dpi=1200)



