

get_ipython().magic('reset -sf')

#code to generat Figure 3 in the main text and SI figure S5-S6

# Imports
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from matplotlib.lines import Line2D
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples


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
data = data[headers2]

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

#kmeans clustering based off the ideal number of clusters
kmeans = KMeans(n_clusters=3, random_state=0)

#need to make sure that the X used here is the same as in classify mutagens
X = StandardScaler().fit_transform(data[headers2])
X = X[:, ~np.isnan(X ).any(axis=0)]
X = X[:, ~np.isinf(X ).any(axis=0)]

principalComponents = PCA(n_components=10).fit_transform(X)
PCA_components = pd.DataFrame(principalComponents)
X = PCA_components.iloc[:,:5]


cluster_labels = kmeans.fit_predict(X)
clusters = kmeans.fit_predict(X)
clust_DF = pd.DataFrame({'cluster':clusters})

# clust_DF2 = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\cluster_data.xlsx", sheet_name = 'cluster_assignments')

filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'cluster_data.xlsx'))

clust_DF2 = pd.read_excel(filename, sheet_name = 'cluster_assignments')



fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)

#larger molecules
ax.scatter(PCA_components[0][clust_DF['cluster'] == 0],
           PCA_components[1][clust_DF['cluster'] == 0], 
           marker ='s',alpha=.8, c = '#004488')  #blue           
           
#smaller molecules
ax.scatter(PCA_components[0][clust_DF['cluster'] == 1],
           PCA_components[1][clust_DF['cluster'] == 1], 
           marker ='v',alpha=.8, c = '#DDAA33')  #  gold

ax.scatter(PCA_components[0][clust_DF['cluster'] == 2],
           PCA_components[1][clust_DF['cluster'] == 2], 
           marker ='o',alpha=.8, c = '#BB5566')   #rose



plt.show()

 
#clean up the axis labels
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.spines['top'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
plt.xlabel('Principle Component 1', fontsize = '14')
plt.ylabel('Principle Component 2',  fontsize = '14')

custom_lines = [Line2D([0], [0], color='#DDAA33', 
                       markersize = 10,  marker = 'v',lw=0, alpha = 0.75),  
    
                Line2D([0], [0], color='#BB5566',linewidth = 1, 
                       markersize = 10,  marker = 'o',lw = 0, alpha = 0.9), 

                Line2D([0], [0], color='#004488', 
                       markersize = 10,  marker = 's',lw=0, alpha = 0.75)] 



ax.legend(custom_lines, ( 'Other Molecules Cluster', 
                         'Large Molecules Cluster', 
                         'Small Molecules Cluster'))

plt.show()

fig.savefig(r'F3_S5_pcaclusters.png', format='png', dpi=1200)


#set up the sillhouette score plot
fig2 = plt.figure()
ax20 = fig2.add_subplot(111)

# Compute the silhouette scores for each sample
silhouette_avg = silhouette_score(X, cluster_labels)
sample_silhouette_values = silhouette_samples(X, cluster_labels)

y_lower = 10
n_clusters = 3

colorMap = ['#004488','#DDAA33','#BB5566' ]
Cluster_Labels = [ 'Small Molecules','Other Molecules', 'Large Molecules']


for i in [1, 0,2]:
    print(i)
    print(len(sample_silhouette_values[clust_DF['cluster'] == i]))
    # Aggregate the silhouette scores for samples belonging to
    # cluster i, and sort them
    ith_cluster_silhouette_values = \
        sample_silhouette_values[cluster_labels == i]

    ith_cluster_silhouette_values.sort()

    size_cluster_i = ith_cluster_silhouette_values.shape[0]
    y_upper = y_lower + size_cluster_i


    #use the same colors as the scatter plot above
    color = cm.nipy_spectral(float(i) / n_clusters)
    
    color = colorMap[i]
    ax20.fill_betweenx(np.arange(y_lower, y_upper),
                      0, ith_cluster_silhouette_values,
                      facecolor=color, edgecolor=color, alpha=0.8)

    # Label the silhouette plots with their cluster numbers at the middle
    ax20.text(0.01, y_lower + 0.5 * size_cluster_i-5, Cluster_Labels[i])
    
    # Compute the new y_lower for next plot
    y_lower = y_upper + 10  # 10 for the 0 samples

#X label
ax20.set_xlabel("Silhouette coefficient values")

# The vertical line for average silhouette score of all the values
ax20.axvline(x=silhouette_avg, color="red", linestyle="--")

# Clear the yaxis labels / ticks
ax20.set_yticks([])  
ax20.set_xticks([-0.1, 0, 0.2, 0.4, 0.6])
plt.grid(None)
plt.show()

fig2.savefig(r'S6_silhouettebladesPCA.png', format='png', dpi=1200)


