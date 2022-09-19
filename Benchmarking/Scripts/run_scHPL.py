import os
from sys import argv
from pathlib import Path
import numpy as np
import pandas as pd
import time as tm
from scHPL import progressive_learning, predict, evaluate
import rpy2.robjects as robjects

def run_scHPL(DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath = "", NumGenes = 0):
    '''
    run scHPL
    Wrapper script to run scHPL on a benchmark dataset with 5-fold cross validation,
    outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
    as row names and gene names as column names.
    LabelsPath : Cell population annotations file path (.csv).
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    OutputDir : Output directory defining the path of the exported file.
    GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
    defining the genes order for each cross validation fold, default is NULL.
    NumGenes : Number of genes used in case of feature selection (integer), default is 0.
    '''
    
    print("read the Rdata file")
    # read the Rdata file
    robjects.r['load'](CV_RDataPath)
    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1
    test_ind = np.array(robjects.r['Test_Idx'])
    train_ind = np.array(robjects.r['Train_Idx'])
    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)
    
    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]
    
    # folder with results
    os.chdir(OutputDir)
    
    train_time = []
    test_time = []
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]
        start = tm.time()
        classifier = 'svm'
        dimred = False
        threshold = 0.25
        tree = progressive_learning.learn_tree([train], [y_train.to_numpy()], classifier = classifier, dimred = dimred, threshold = threshold)
        train_time.append(tm.time()-start)
	    # predict
        start = tm.time()
        ypred = predict.predict_labels(test, tree)
        test_time.append(tm.time()-start)
        truelab.extend(y_test.values)
        pred.extend(ypred)
        print(pred)
    
            
    true = pd.DataFrame(truelab)
    prediction = pd.DataFrame(pred)
    tr_time = pd.DataFrame(train_time)
    ts_time = pd.DataFrame(test_time)
    
    os.chdir(OutputDir)
    true.to_csv("scHPL_true.csv",index = False)
    prediction.to_csv("scHPL_pred.csv",index = False)
    tr_time.to_csv("scHPL_training_time.csv",index = False)
    ts_time.to_csv("scHPL_test_time.csv",index = False)



        
        
run_scHPL(argv[1], argv[2], argv[3], argv[4])
        
        
        






