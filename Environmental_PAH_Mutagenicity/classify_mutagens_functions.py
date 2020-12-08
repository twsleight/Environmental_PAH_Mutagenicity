#Gaussian File Reader and Structural Descriptor Anaylsis
#Updated on 5.21.2020
#Trevor Sleight and Caitlin Sexton
#****************************************************

# import os
import pandas as pd
import numpy as np


from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix 

import time
from sklearn import metrics
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant

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
        
        
        #diagnostic variable, eliminate later        
        # columnNames = list(finaldata.columns)   
        # [a] = logmodel.coef_
        # rankedColumns_Raw = pd.DataFrame(data = {'Coef':a, 'Name':columnNames})
         
        #***********************************                
        #need to change this to predict proba in order to be able to 
        #manipulate the threshold for better recall
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
def removehighestVIF(finaldata, VIF_threshhold):
    
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
    
    worstDesc = vifData.idxmax()
    #remove this descriptor from finaldata
    finaldata = finaldata.drop(columns = [worstDesc])
    print('dropping ', worstDesc)

    
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
        finaldata = RemoveHighestVIF(finaldata, VIF_threshhold)
    #if you get to here, something has gone seriously wrong. 

    return finaldata
####################################


