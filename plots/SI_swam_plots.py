# from IPython import get_ipython;   
# get_ipython().magic('reset -sf')

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



# 0 is the small cluster of 240 molecules
# 2 is the large molecules cluster of 315 molecules

#use to select data from just one cluster
data = data.loc[clust_DF['cluster']==2]

#set up the plotting parameters

# plot_Data = pd.DataFrame(data = {'HOMO':data['HOMO'],
#                                  'HOMO' : data['HOMO']*27.2114,                                 
#                                  'Mutagenicity':data['result'] } )

plot_Data = pd.DataFrame(data = {'SpMin5_Bhi':data['SpMin5_Bhi'],
                                  'SpMin5_Bhi' : data['SpMin5_Bhi'],                                 
                                  'Mutagenicity':data['result'] } )

# plot_Data = pd.DataFrame(data = {'RDF55s':data['RDF55s'],
#                                   'RDF55s' : data['RDF55s']*27.2114,                                 
#                                   'Mutagenicity':data['result'] } )



plot_Data.replace({1: 'Positive' , 0:'Negative' }, inplace = True)


fig = plt.figure(figsize = (10,8))
ax1 = fig.add_subplot(1, 1, 1)

sns.swarmplot(x="Mutagenicity", y="SpMin5_Bhi", data=plot_Data, ax = ax1)

sns.set_theme(style="whitegrid")
sns.despine()

# plt.ylabel("HOMO-LUMO Gap (eV)", fontsize = 18, labelpad = 15)

plt.ylabel("SpMin5_Bhi", fontsize = 18, labelpad = 8)
plt.xlabel("Mutagenicity", fontsize = 18, labelpad = 15)


plt.yticks(fontsize=15) 
plt.xticks(fontsize=15) 


# fig.savefig(r'SpMin5_Bhi_swam_all.jpg', format='jpg', dpi=1200)