from IPython import get_ipython;   
get_ipython().magic('reset -sf')

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os


parent = os.path.join(os.path.abspath(__file__), os.pardir)

filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'mutagenicity_data.xlsx'))

data = pd.read_excel(filename, sheet_name = 'Sheet1')


filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'cluster_data.xlsx'))

clust_DF = pd.read_excel(filename, sheet_name = 'cluster_assignments')



# 0 is the small cluster of 215 molecules
# 2 is the large molecules cluster of 301 molecules

#use to select data from just one cluster
data = data.loc[clust_DF['cluster']==0]

#set up the plotting parameters

plot_Data = pd.DataFrame(data = {'RDF55s':data['RDF55s'],
                                 'HOMO-LUMO Gap (eV)' : data['HLgap'],                                 
                                 'Mutagenicity':data['result'] } )

plot_Data.replace({1: 'Positive' , 0:'Negative' }, inplace = True)


fig = plt.figure(figsize = (10,8))
ax1 = fig.add_subplot(1, 1, 1)

sns.swarmplot(x="Mutagenicity", y="RDF55s", data=plot_Data, ax = ax1)

sns.set_theme(style="whitegrid")
sns.despine()

# plt.ylabel("HOMO-LUMO Gap (eV)", fontsize = 18, labelpad = 15)

plt.ylabel("RDF55s", fontsize = 18, labelpad = 8)
plt.xlabel("Mutagenicity", fontsize = 18, labelpad = 15)


plt.yticks(fontsize=15) 
plt.xticks(fontsize=15) 


# fig.savefig(r'RDF55s_swam_all.jpg', format='jpg', dpi=1200)