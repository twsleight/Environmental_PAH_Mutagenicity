#Gaussian File Reader and Structural Descriptor Anaylsis
#Updated on 5.21.2020
#Trevor Sleight and Caitlin Sexton
#****************************************************

from IPython import get_ipython;   
get_ipython().magic('reset -sf')

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import gridspec

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix 

import time
from sklearn import metrics
from sklearn.decomposition import PCA

from sklearn.metrics import auc
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant
import copy

import xlsxwriter
from collections import Counter

# Run classifier with cross-validation and plot ROC curves
#supporting functions. 
def regress_kfold(finaldata, y, n_splits=3):
    '''impletments a k-fold cross validation and returns a confusion matrix

    Args:finaldata, y, n_splits
    Returns:conf
    '''    

    tprs = []
    aucs = []
    threshs = []
    confs = np.array([[0,  0],[ 0, 0]])
    mean_fpr = np.linspace(0, 1, 50)
         
    cv = StratifiedKFold(n_splits=n_splits)

    for i, (train, test) in enumerate(cv.split(finaldata, y)):
        
        X_train = finaldata.iloc[train,:]
        X_test = finaldata.iloc[test,:]    
        y_train = y[train]
        y_test = y[test]  
    
        #get a new random seed from the clock
        cValue = int(time.perf_counter()*1000)
        
        #This needs to be inside the loop because we need to reinitalize it 
        #with the clock value each time.         
        logmodel = LogisticRegression(multi_class = 'ovr',
                                      solver = 'liblinear', 
                                      random_state=cValue, 
                                      max_iter = 1000000, 
                                      class_weight={1:5, 0:1})
        logmodel.fit(X_train, y_train)
        

        y_test_probabilities = logmodel.predict_proba(X_test)
        predictions = (y_test_probabilities[:,1] > .5).astype(int)

        confs += confusion_matrix(y_test, predictions)
        
        y_pred_proba = logmodel.predict_proba(X_test)[::,1]
        
        fpr, tpr, thresholds = metrics.roc_curve(y_test,  y_pred_proba, pos_label = 1)
        auc1 = metrics.roc_auc_score(y_test, y_pred_proba)
    
        #need to use the np interp because the exact size of the data is 
        #going to vary with the CV splits
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_threshs = np.interp(mean_fpr, fpr, thresholds)
               
        tprs.append(interp_tpr)               
        threshs.append(interp_threshs)
        aucs.append(auc1)
       
    return (confs, tprs, aucs, threshs)

#############################################
def RemoveHighestVIF(finaldata, VIF_threshhold, pref_descriptors):
    
    '''recursively removes correlated descriptors based on VIF
    until all descriptors are below the specified VIF threshold.
    
    Optional: prioritized descriptors based on a list. 

    Args:finaldata, VIF_threshold
    
    Returns:finaldata
    '''   

    #find the highest VIF
    #***********************
    X = add_constant(finaldata)    
    vifData = pd.Series([variance_inflation_factor(X.values, i) 
               for i in range(X.shape[1])], index=X.columns)
       
    #remove const and then find the worst descriptor
    vifData.drop('const', inplace = True)
    

    #look at the two worst descriptors. 
    #these will be correlated. if they're within 0.5 VIF of each other, drop
    #the less preferred descriptor. If there's a greater spread than that
    #drop the highest VIF value
      
    ####################################
    vifData.sort_values(inplace = True, ascending = False)
    spread_vif = vifData[0]-vifData[1]
    
    print('\nthe VIF spread is ', spread_vif)
    # print(vifData.index[0], vifData.index[1])
 
    #if the spread of the VIF is less that
    if spread_vif < .5 or np.isnan(spread_vif) :
        #choose the preferred descriptor
        for i, row in pref_descriptors.iterrows():
            if vifData.index[0] in row['Descriptor']:
                rank0Desc = i

            #needs to be a seperate if statment 
            #since they could be the same category
            if vifData.index[1] in row['Descriptor']:
                rank1Desc = i
        
        #lower numbers are more preferred
        if rank0Desc > rank1Desc:
            print(' preference selected, keeping ',vifData.index[1] , ' dropping ',vifData.index[0]  )
            worstDesc = vifData.index[0]
        elif rank1Desc > rank0Desc:
            print(' preference selected, keeping ',vifData.index[0], ' dropping ',vifData.index[1]   )            
            worstDesc = vifData.index[1]
        elif rank0Desc == rank0Desc:
            worstDesc = vifData.idxmax()
            print('dropping ', worstDesc,' keeping ', vifData.index[1])        
        else:
            raise "Error"
    else:
        worstDesc = vifData.idxmax()
        print('dropping ', worstDesc,' keeping ', vifData.index[1])        
    ####################################

    #remove this descriptor from finaldata
    finaldata = finaldata.drop(columns = [worstDesc])

    #get a new VID dataframe on the updated finaldata dataframe
    X = add_constant(finaldata)    
    vifData = pd.Series([variance_inflation_factor(X.values, i) 
               for i in range(X.shape[1])], 
              index=X.columns)
    
    #if we're now under the threshold we're done
    if all(vifData < VIF_threshhold):
        print('should be returning')
        #now we're done,return the function 
        return finaldata
    #otherwise go around again
    else:
        print('got here')
        finaldata = RemoveHighestVIF(finaldata, VIF_threshhold, pref_descriptors)
    #if you get to here, something has gone seriously wrong. 

    return finaldata
####################################

#changes working directory
os.chdir(r'C:\ResearchWorkingDirectory\finalData_QSAR\rawData')

# data = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\mutagenicity_data.xlsx", sheet_name = 'Sheet1')

data = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\small_mutagenicity_data.xlsx", sheet_name = 'Sheet1')

#create a copy that won't be modified by any processing later
headersx = list(data.columns)
dataOrig = copy.deepcopy(data[headersx])

headers2 = list(data.columns)
headers2.remove('filename')
headers2.remove('IP')
headers2.remove('EA')
headers2.remove('Neutral')
headers2.remove('SMILES')
headers2.remove('result')
headers2.remove('CAS')
headers2.remove('name')
# headers2.remove('Unnamed: 0')

dropnondata = []
for h in headers2: 
    if sum(data[h]) == 0:
        dropnondata.append(h)
    # if nan in np.array(data[h]):
    #     print('nan found')

headers2 = list(set(headers2) - set(dropnondata))

#do the clusting without considering the truth data
data = data[headers2]
#PCA clustering
X = StandardScaler().fit_transform(data[headers2])
X = X[:, ~np.isnan(X ).any(axis=0)]
X = X[:, ~np.isinf(X ).any(axis=0)]

# len(headers2)


principalComponents = PCA(n_components=10).fit_transform(X)
PCA_components = pd.DataFrame(principalComponents)
X = PCA_components.iloc[:,:5]

#initialized the random state so we get consistent results
kmeans = KMeans(n_clusters=3, random_state=0)
clusters = kmeans.fit_predict(X)
clust_DF = pd.DataFrame({'cluster':clusters, 'SMILES': dataOrig['SMILES']}) 


clust_DF = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\cluster_data.xlsx", sheet_name = 'cluster_assignments')

#use for data accross the full equation
data = dataOrig

#use to select data from just one cluster
# data = dataOrig.loc[clust_DF['cluster']==0]

# file 1 goes with cluster 0, small molecules
#cluster 1 is other molecules
#clsuter 2 is larger molecules file 0 goes with it. 

#the _1 file goes with the smaller molecules cluster, 210 molecules in it. 
bestData = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\test_DescSelect_all.xlsx", sheet_name = 'Descriptors')


#************************************
#This processes the output from the RFE analysis

allselectDesc = []

#for the example dataset - finalIDS 50
for i in range(30, 901,90):


#for full dataset - finalIDS 72
# for i in range(46, 901,90):

    #larger molecules cluster - finalIDS 45
# for i in range(23, 901,90):
    
    #small molecules - finalIDs 32
# for i in range(16, 901,90):
      
    curDesc = list(bestData[i].dropna()) 

    
    allselectDesc.extend(list(set(curDesc)))

    print(len(set(allselectDesc)))
                            
finalIDS = list(set(allselectDesc))  #['HLgap']


#305 is cluster 0, 40 is cluster 2, 210 is cluster 1

#get the descriptors that were found 3 times or more
a = Counter(allselectDesc)
b = pd.DataFrame(data = {'values':a.values(), 'descriptors':a.keys() } )
# c = b
c = b.loc[b['values'] >= 3]

finalIDS = list(c['descriptors'])
finalIDS = list(set(finalIDS).intersection(set(headers2)))

finaldata = data[finalIDS]

# Classification and ROC analysis section
#randomize the order of data
tprs = []
aucs = []
threshs = []
confs = np.array([[0,  0],[ 0, 0]])
mean_fpr = np.linspace(0, 1, 50)
    

#only use the columns designated by finalIDS from here on
data[finalIDS] = StandardScaler().fit_transform(data[finalIDS] )

for n in range(0,10):
    #shuffle the order of the data rows at the n level
    data2 = data.sample(frac=1).reset_index(drop=True)
    
    finaldata = data2[finalIDS]
    y = data2['result']

    confs_results = regress_kfold(finaldata, y, n_splits=3)
    
    confs = np.add(confs,confs_results[0])
    tprs.extend(confs_results[1])   
    aucs.extend(confs_results[2])
    threshs.extend(confs_results[3])  
             

#plotting section
#***************************************************
tn = confs[0,0]
tp = confs[1,1]
fn = confs[1,0]
fp = confs[0,1]

allTP = tp + fn
allTN = tn + fp
allvals = allTP + allTN

acc = (tp+tn) / (tp + tn + fn + fp)

print(confs)

print('accuracy is ', acc) 
print('Precision is ', tp/(tp+fp))
prec = tp/(tp+fp)
print('Recall/Sensitivity ', tp/(tp+fn))  
rec =  tp/(tp+fn)
print('Specificity is ', tn/(tn+fp))
spec = tn/(tn+fp)
print('F1 is ', 2*(prec*rec)/(prec+rec))

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0

mean_thresh = np.mean(threshs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)




#ROC Plot************************************

# coefficient column on each side
fig = plt.figure(figsize=(13, 8)) 
gs = gridspec.GridSpec(1,6, width_ratios=[.2,.2, 1.95,0.05, 0.2, 0.15]) 
ax0 = plt.subplot(gs[2])
#******************************


ax0.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
        label='Chance', alpha=.8)

ax0.plot(mean_fpr, mean_tpr, color='b',
        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
        lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax0.fill_between(mean_fpr, tprs_lower, tprs_upper,color='blue',alpha=.2, label=r'$\pm$ 1 std. dev.')

ax0.plot(mean_fpr,tprs_lower, color='k',linestyle = '--',lw =1, alpha=.8)
ax0.plot(mean_fpr,tprs_upper, color='k',linestyle = '--',lw =1, alpha=.8)

ax0.set_ylabel('Specificity', fontsize = 16, labelpad = 10)
ax0.set_xlabel('Sensitivity/Recall', fontsize = 16, labelpad = 10)
ax0.tick_params(axis='both', labelsize=14 )
ax0.legend(loc="lower left", fontsize = 12)
ax0.grid(False)


#need to add text for classifier stuff
F1=  2*(prec*rec)/(prec+rec)

textstr = '\n'.join((
    r'$\mathrm{Accuracy}=%.2f$' % (acc, ),
    r'$\mathrm{Precision}=%.2f$' % (prec, ),
    r'$\mathrm{Sensitivity/Recall}=%.2f$' % (rec, ),
    r'$\mathrm{Specificity}=%.2f$' % (spec, ),
    r'$\mathrm{F1}=%.2f$' % (F1, )))

props = dict(boxstyle='round', facecolor='white', alpha=0.5)

# place a text box in upper left in axes coords
ax0.text(0.98,0.4, textstr, transform=ax0.transAxes, fontsize=14,
        verticalalignment='center',horizontalalignment = 'right', bbox=props)

#regress over the entire dataset to get the final regression coefficients. 

#start from 0 so this is consistent
logmodel = LogisticRegression(multi_class = 'ovr',
                              solver = 'liblinear', 
                              random_state=0, 
                              max_iter = 1000000, 
                              class_weight={1:5, 0:1})
logmodel.fit(finaldata, y)

#get the data for the residuals
proba_preds = logmodel.predict_proba(finaldata)[::,1]
resids = y.subtract(proba_preds)


columnNames = list(finaldata.columns)   
[a] = logmodel.coef_
rankedColumns_Raw = pd.DataFrame(data = {'Coef':a, 'Name':columnNames})
# rankedColumns_Raw.replace(to_replace = abbrev_ref,inplace = True)

dfx = rankedColumns_Raw.set_index(keys = 'Name', drop = True).sort_values(by = 'Coef', ascending = False)


#setup for a 2 sized heatmap
import matplotlib
matrix = np.array([0,1.5])
boundaries = [value for value in matrix.flatten().tolist()]
list.sort(boundaries)
colors = ["#0000ff", "#FFFF00"]
norm = matplotlib.colors.BoundaryNorm(boundaries=boundaries + [boundaries[-1]], ncolors=256)

cmapR = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

matrix = np.array([-1.5,0])
boundaries = [value for value in matrix.flatten().tolist()]
list.sort(boundaries)
colors = ["#FFFF00", "#ff0000"]
norm = matplotlib.colors.BoundaryNorm(boundaries=boundaries + [boundaries[-1]], ncolors=256)

cmapL = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)


#*************column on each side***********8

dfxNeg = dfx.loc[dfx['Coef'] <0]
dfxPos = dfx.loc[dfx['Coef'] >0]
maxVal = max(dfx['Coef'])
minVal = min(dfx['Coef'])
spreadVal = min(dfxPos['Coef']) - max(dfxNeg['Coef'])


ax1 = plt.subplot(gs[0])
sns.heatmap( dfxPos,ax = ax1, vmin = minVal+spreadVal, vmax = maxVal, cmap="jet", linewidths=0.5, 
            annot=True, cbar = False, fmt = ".2f", 
            annot_kws={"size":"12","fontweight":"bold"}, yticklabels=True)

ax1.tick_params(labelleft = True, labelright = False, rotation=0)
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.set_anchor('W')
ax1.set_title('Positive Coefficients', fontsize = 12, pad = 15)
ax1.tick_params(left=False, bottom=False, labelsize = 12)


ax2 = plt.subplot(gs[4])
sns.heatmap( dfxNeg,ax = ax2,vmin = minVal, vmax = maxVal-spreadVal, cmap="jet", linewidths=0.5, 
            annot=True, cbar = False, fmt = ".2f", 
            annot_kws={"size":"12","fontweight":"bold"}, yticklabels=True)


ax2.tick_params(labelleft = False, labelright = True, rotation=0)
ax2.set_xlabel("")
ax2.set_ylabel("")
ax2.set_anchor('W')
ax2.set_title('Negative Coefficients', fontsize = 12, pad = 15)
ax2.tick_params(left=False, bottom=False, labelsize = 12)

plt.show()

#Save a plot at this point if you'd like. 
# fig.savefig(r'C:\Users\twsle\OneDrive\Documents\Second Paper\images\Final\ROC_all_vif.jpg', format='jpg', dpi=1200)


#VIF section
###################################
X = add_constant(finaldata)
vifData = pd.Series([variance_inflation_factor(X.values, i) 
               for i in range(X.shape[1])], index=X.columns)

pref_descriptors = pd.read_excel(r"C:\ResearchWorkingDirectory\Environmental_PAH_Mutagenicity\Final_Data\DescriptorPriority.xlsx", sheet_name = 'Sheet1')

dataX = finaldata
dataX = RemoveHighestVIF(finaldata,5,pref_descriptors)
finalColumns = dataX.columns

initalColumns = finaldata.columns
correlations = finaldata[initalColumns].corr(method='spearman')

X = add_constant(dataX)
vifData2 = pd.Series([variance_inflation_factor(X.values, i) 
               for i in range(X.shape[1])], index=X.columns)
###################################


# plot correlation matrix - only do this before the VIF removes
#******************************************
fig2 = plt.figure(figsize = (8,8))
plt.rcParams['figure.constrained_layout.use'] = False

ax = fig2.add_subplot(111)
cax = ax.matshow(correlations, vmin=-1, vmax=1, cmap = 'jet')

bar = fig2.colorbar(cax,shrink =0.9)

fig2.canvas.flush_events() #else bar.ax.get_yticklabels() is not yet updated
bar.ax.set_yticklabels(labels=bar.ax.get_yticklabels(), fontsize=12)  #, weight='bold',

ticks = np.arange(0,len(initalColumns),1)
ax.set_xticks(ticks)
ax.set_yticks(ticks)

ax.set_xticklabels(list(initalColumns), fontsize = 10, rotation=-90, 
                   ha = 'left', rotation_mode="anchor", fontweight = "bold")

ax.xaxis.set_ticks_position('bottom')

ax.set_yticklabels(list(initalColumns), fontsize = 10, fontweight = "bold")
ax.grid(False)

fig2.subplots_adjust(bottom=0.2, left = 0.2, top =0.9, right = 1)
plt.show()
#******************************************************************

#save off the correlation image
# fig2.savefig(r'C:\Users\twsle\OneDrive\Documents\Second Paper\images\Final\Correlation_large.jpg', format='jpg', dpi=1200)
#******************************************************************

#redo the analysis for the ROC plot from above with  VIF reduced data

finaldata = dataX
finalIDS = list(dataX.columns)   

   
finaldata = data[finalIDS]

# Classification and ROC analysis section

#randomize the order of data
tprs = []
aucs = []
threshs = []
confs = np.array([[0,  0],[ 0, 0]])
mean_fpr = np.linspace(0, 1, 50)
    

#only use the columns designated by finalIDS from here on
data[finalIDS] = StandardScaler().fit_transform(data[finalIDS] )

for n in range(0,10):
    #shuffle the order of the data rows at the n level
    data2 = data.sample(frac=1).reset_index(drop=True)
    
    finaldata = data2[finalIDS]
    y = data2['result']

    confs_results = regress_kfold(finaldata, y, n_splits=3)
    
    confs = np.add(confs,confs_results[0])
    tprs.extend(confs_results[1])   
    aucs.extend(confs_results[2])
    threshs.extend(confs_results[3])  
             

#plotting section
#***************************************************
tn = confs[0,0]
tp = confs[1,1]
fn = confs[1,0]
fp = confs[0,1]

allTP = tp + fn
allTN = tn + fp
allvals = allTP + allTN

acc = (tp+tn) / (tp + tn + fn + fp)

print(confs)

print('accuracy is ', acc) 
print('Precision is ', tp/(tp+fp))
prec = tp/(tp+fp)
print('Recall/Sensitivity ', tp/(tp+fn))  
rec =  tp/(tp+fn)
print('Specificity is ', tn/(tn+fp))
spec = tn/(tn+fp)
print('F1 is ', 2*(prec*rec)/(prec+rec))

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0

mean_thresh = np.mean(threshs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)




#ROC Plot************************************

# coefficient column on each side
fig = plt.figure(figsize=(13, 8)) 
gs = gridspec.GridSpec(1,6, width_ratios=[.2,.2, 1.95,0.05, 0.2, 0.15]) 
ax0 = plt.subplot(gs[2])
#******************************


ax0.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
        label='Chance', alpha=.8)

ax0.plot(mean_fpr, mean_tpr, color='b',
        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
        lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax0.fill_between(mean_fpr, tprs_lower, tprs_upper,color='blue',alpha=.2, label=r'$\pm$ 1 std. dev.')

ax0.plot(mean_fpr,tprs_lower, color='k',linestyle = '--',lw =1, alpha=.8)
ax0.plot(mean_fpr,tprs_upper, color='k',linestyle = '--',lw =1, alpha=.8)

ax0.set_ylabel('Specificity', fontsize = 16, labelpad = 10)
ax0.set_xlabel('Sensitivity/Recall', fontsize = 16, labelpad = 10)
ax0.tick_params(axis='both', labelsize=14 )
ax0.legend(loc="lower left", fontsize = 12)
ax0.grid(False)


#need to add text for classifier stuff
F1=  2*(prec*rec)/(prec+rec)

textstr = '\n'.join((
    r'$\mathrm{Accuracy}=%.2f$' % (acc, ),
    r'$\mathrm{Precision}=%.2f$' % (prec, ),
    r'$\mathrm{Sensitivity/Recall}=%.2f$' % (rec, ),
    r'$\mathrm{Specificity}=%.2f$' % (spec, ),
    r'$\mathrm{F1}=%.2f$' % (F1, )))

props = dict(boxstyle='round', facecolor='white', alpha=0.5)

# place a text box in upper left in axes coords
ax0.text(0.98,0.4, textstr, transform=ax0.transAxes, fontsize=14,
        verticalalignment='center',horizontalalignment = 'right', bbox=props)

#regress over the entire dataset to get the final regression coefficients. 

#start from 0 so this is consistent
logmodel = LogisticRegression(multi_class = 'ovr',
                              solver = 'liblinear', 
                              random_state=0, 
                              max_iter = 1000000, 
                              class_weight={1:5, 0:1})
logmodel.fit(finaldata, y)

#get the data for the residuals
proba_preds = logmodel.predict_proba(finaldata)[::,1]
resids = y.subtract(proba_preds)


columnNames = list(finaldata.columns)   
[a] = logmodel.coef_
rankedColumns_Raw = pd.DataFrame(data = {'Coef':a, 'Name':columnNames})
# rankedColumns_Raw.replace(to_replace = abbrev_ref,inplace = True)

dfx = rankedColumns_Raw.set_index(keys = 'Name', drop = True).sort_values(by = 'Coef', ascending = False)


#setup for a 2 sized heatmap
import matplotlib
matrix = np.array([0,1.5])
boundaries = [value for value in matrix.flatten().tolist()]
list.sort(boundaries)
colors = ["#0000ff", "#FFFF00"]
norm = matplotlib.colors.BoundaryNorm(boundaries=boundaries + [boundaries[-1]], ncolors=256)

cmapR = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

matrix = np.array([-1.5,0])
boundaries = [value for value in matrix.flatten().tolist()]
list.sort(boundaries)
colors = ["#FFFF00", "#ff0000"]
norm = matplotlib.colors.BoundaryNorm(boundaries=boundaries + [boundaries[-1]], ncolors=256)

cmapL = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)


#*************column on each side***********8

dfxNeg = dfx.loc[dfx['Coef'] <0]
dfxPos = dfx.loc[dfx['Coef'] >0]
maxVal = max(dfx['Coef'])
minVal = min(dfx['Coef'])
spreadVal = min(dfxPos['Coef']) - max(dfxNeg['Coef'])


ax1 = plt.subplot(gs[0])
sns.heatmap( dfxPos,ax = ax1, vmin = minVal+spreadVal, vmax = maxVal, cmap="jet", linewidths=0.5, 
            annot=True, cbar = False, fmt = ".2f", 
            annot_kws={"size":"12","fontweight":"bold"}, yticklabels=True)

ax1.tick_params(labelleft = True, labelright = False, rotation=0)
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.set_anchor('W')
ax1.set_title('Positive Coefficients', fontsize = 12, pad = 15)
ax1.tick_params(left=False, bottom=False, labelsize = 12)


ax2 = plt.subplot(gs[4])
sns.heatmap( dfxNeg,ax = ax2,vmin = minVal, vmax = maxVal-spreadVal, cmap="jet", linewidths=0.5, 
            annot=True, cbar = False, fmt = ".2f", 
            annot_kws={"size":"12","fontweight":"bold"}, yticklabels=True)


ax2.tick_params(labelleft = False, labelright = True, rotation=0)
ax2.set_xlabel("")
ax2.set_ylabel("")
ax2.set_anchor('W')
ax2.set_title('Negative Coefficients', fontsize = 12, pad = 15)
ax2.tick_params(left=False, bottom=False, labelsize = 12)

plt.show()

#this is the final plot
# fig.savefig(r'C:\Users\twsle\OneDrive\Documents\Second Paper\images\Final\ROC_all_vif.jpg', format='jpg', dpi=1200)


#******************************************************************
# #save the fitted model for use later
# from joblib import dump
# dump(logmodel, 'smaller.joblib') 

# dump(kmeans, 'kmeans_clusters.joblib') 


#also save off the columns that should be used with the model
# df = pd.DataFrame(columnNames)
# writer = pd.ExcelWriter(r'C:\ResearchWorkingDirectory\Dissertation_Specific_Code\SmallerColumnNames.xlsx', engine='xlsxwriter')
# df.to_excel(writer, sheet_name='Sheet1', index=False)
# writer.save()


