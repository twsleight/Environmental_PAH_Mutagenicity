
import pandas as pd
import os
import time

import copy
import random

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import RFE
import matplotlib.pyplot as plt 

from sklearn.cluster import KMeans
from classify_mutagens import regress_kfold


#script for selecting the descriptors. Outputs which descriptors were selected
#by recursive feature elimination and the F1 scores that correspone
#This is a time consuming process and usually takes several hours. It 
#is recommended that this be done on a decidated system such as a cluster

plt.rc("font", size=14)
start_time = time.time()

parent = os.path.join(os.path.abspath(__file__), os.pardir)
filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'mutagenicity_data.xlsx'))

data = pd.read_excel(filename, sheet_name = 'Sheet1')

headersx = list(data.columns)
dataOrig = copy.deepcopy(data[headersx])

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

dropnondata = []
for h in headers2: 
    if sum(data[h]) == 0:
        dropnondata.append(h)

headers2 = list(set(headers2) - set(dropnondata))

X = data[headers2]
X = StandardScaler().fit_transform(X)

#check the array for infs or NaNs

#need to adjust the clustering to use PCA instead of raw data and see if it changes anything. 
pca = PCA(n_components=10)
principalComponents = pca.fit_transform(X)
PCA_components = pd.DataFrame(principalComponents)

X = PCA_components.iloc[:,:5]

kmeans = KMeans(n_clusters=3, random_state=0)
clusters = kmeans.fit_predict(X)
clust_DF = pd.DataFrame({'cluster':clusters,'result': data['result'], 'SMILES': dataOrig['SMILES']})

data = data.loc[clust_DF['cluster']==0]

print(len(data['SMILES']))

headers3 = copy.deepcopy(headers2)
headers3.append('result')

#scale the descriptors, but not the truth values
y = list(data['result'])
data[headers2] = StandardScaler().fit_transform(data[headers2]) 
data['result'] = y
#from this point we don't touch the data again. leave it pristine


logmodel = LogisticRegression(multi_class = 'ovr',solver = 'liblinear', random_state=1, max_iter = 1000000)

allF1Scores = pd.DataFrame()
selectedDesc = pd.DataFrame()

for i in range(0,10):
    print(i)
    #use the same clock value for each of the possible numbers of features so that these are comparable. 
    cValue = int(str(time.perf_counter()*1000)[:4])    
    
    
    #shuffule the order of the data and the headers, need to make sure to shuffle the answers also
    data2 = data.sample(frac=1).reset_index(drop=True)    
    random.shuffle(headers2)
       
    #split the truth data and the deccriptors
    y = data2['result']   
    finaldata = data2[headers2]

    F1Scores = pd.DataFrame(columns = ['F1'])  

    #custom recursive feature elimination      
    for nfeat in range(10, 100, 1):
        print('the number of features is ', nfeat)
        
        #RFE works over the entire dataset
        selector = RFE(estimator = logmodel, n_features_to_select = nfeat, step = 10)
        selector = selector.fit(finaldata, y)
        rfe_fits = selector.ranking_
        
        columnNames = finaldata.columns
        rankedColumnns_Raw = pd.DataFrame(data = {'Rank':selector.ranking_, 'Name':columnNames})
        
                
        #data to use
        evalData = list(rankedColumnns_Raw[rankedColumnns_Raw['Rank']==1]['Name'] )
            
        tempSelDesc = rankedColumnns_Raw.loc[rankedColumnns_Raw['Rank']==1].reset_index(drop = True).drop(labels = ['Rank'], axis = 1)
        tempSelDesc = tempSelDesc.rename(columns={'Name': nfeat})
        
        #add the data from this iteration to the existing data
        selectedDesc = pd.concat([selectedDesc,tempSelDesc], ignore_index = True, sort = False, axis = 1)
        
        
     
        print('the number of descriptors is ',len(evalData))
        
        #this is the data that we're using
        selectdata = data2[evalData]

        # X_train, X_test, y_train, y_test = train_test_split(finaldata, y, test_size=0.33, random_state=cValue, stratify = y)

#         logmodel = LogisticRegression(multi_class = 'ovr',solver = 'liblinear', random_state=1, max_iter = 1000000)
#         logmodel.fit(X_train, y_train)
        
#         print(sum(y_test))
        
#         y_test_probabilities = logmodel.predict_proba(X_test)
#         predictions = (y_test_probabilities[:,1] > 0.5).astype(int)
# #       predictions = logmodel.predict(X_test, cutoff = 0.4)
        
#         print(confusion_matrix(y_test, predictions))
# #       print(accuracy_score(y_test, predictions))
#         conf = confusion_matrix(y_test, predictions)  
      
        
      
        confs_results = regress_kfold(selectdata, y, n_splits=3)        
        conf = confs_results[0]

        tn = conf[0,0]
        tp = conf[1,1]
        fn = conf[1,0]
        fp = conf[0,1]
        
        acc = (tp+tn) / (tp + tn + fn + fp)                       
        prec = tp/(tp+fp)
        rec =  tp/(tp+fn)   
                  
#            print('accuracy is ', acc) 
#            print('Precision is ', tp/(tp+fp))                                   
#            print('Recall/Sensitivity ', tp/(tp+fn))  
#            print('Specificity is ', tn/(tn+fp))          
#            print('F1 is ', 2*(prec*rec)/(prec+rec))
        
        tempSlice =  pd.DataFrame(data = { 'F1': 2*(prec*rec)/(prec+rec)}, index =[nfeat] )
        
        F1Scores = F1Scores.append(tempSlice)
        #now I need a method of checking how well the selected features evaluted the data
            
    #lets join the f1 dataframes on the nfeat index
    allF1Scores = pd.concat([allF1Scores, F1Scores,], axis=1)


#put the F1 scores in one dataframe and the the descriptors in another, put these in two spreadsheet tabs
#allF1Scores.to_excel("test_DescSelect.xlsx", sheet_name='F1')  
#selectedDesc.to_excel("test_DescSelect.xlsx", sheet_name='Descriptors')  

print("--- %s seconds ---" % (time.time() - start_time))


#write off the data - adjust the files names as appropriate

# writer = pd.ExcelWriter('test_DescSelect_all.xlsx')
# selectedDesc.to_excel(writer,'Descriptors')
# allF1Scores.to_excel(writer,'F1')
# writer.save()


# writer = pd.ExcelWriter('cluster_data.xlsx')
# clust_DF.to_excel(writer,'cluster_assignments')
# writer.save()



