import os 
import numpy as np
import pandas as pd
import time as tm
from sklearn.svm import LinearSVC
from sys import argv
from pathlib import Path
from sklearn.calibration import CalibratedClassifierCV


def run_SVMrej(RefPath, LabelsPath, TestPath, OutputDir, rejection):
	'''
	run baseline classifier: SVM
	Wrapper script to run an SVM classifier with a linear kernel 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath : Ref Data file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv).
	TestPath : Test dataset path
	OutputDir : Output directory defining the path of the exported file.
	'''
	os.system("/project/kleinman/hussein.lakkis/from_hydra/2021_07_15-Pipeline/Scripts/MRMR. R {0} {1} 2000 {2}".format(RefPath, LabelsPath, OutputDir))
	# read the data
	data = pd.read_csv(RefPath,index_col=0,sep=',')
	labels = pd.read_csv(LabelsPath, index_col = 0)
	test = pd.read_csv(TestPath,index_col=0,sep=',')
	data.columns = map(str.upper, data.columns)
	test.columns = map(str.upper, test.columns)
	common_genes = np.intersect1d(data.columns, test.columns).tolist()
	test = test[common_genes]
	data = data[common_genes]
	# this removes dups
	data = data.groupby(level=0, axis=1).first()
	test = test.groupby(level=0, axis=1).first()
	# Add assertion statement
	# get the Y vector (targets)
	y_train = np.array(labels['label'])
	# folder with results
	os.chdir(OutputDir)
	# normalize data
	data = np.log1p(data)
	clf = LinearSVC()
	Classifier = CalibratedClassifierCV(clf)

	tr_time=[]
	ts_time=[]
	start=tm.time()
	Classifier.fit(data, y_train)
	tr_time.append(tm.time()-start)

	start=tm.time()
	predicted = Classifier.predict(test)
	# predict the probabilities 
	prob = Classifier.predict_proba(test)
	ts_time.append(tm.time()-start)
	probabilities = np.max(prob, axis = 1)
	# create a dataframe with the output
	if reject:
		unlabeled = np.where(probabilities < threshold)
		predicted[unlabeled] = 'Unknown'
		pred = pd.DataFrame({'SVM_Predictions_Rej' : predicted,'SVM_Predictions':Classifier.predict(test), 'SVM_Probability_score' : probabilities})
		pred.index = test.index
	else:
		pred = pd.DataFrame({'SVM_Predictions' : predicted, 'SVM_Probability_score' : probabilities})


	ts_time.append(tm.time()-start)

	tr_time = pd.DataFrame(tr_time)
	ts_time = pd.DataFrame(ts_time)

	pred.to_csv("SVM_reject_pred.csv", index = True)
	tr_time.to_csv("SVM_reject_training_time.csv", index = False)
	ts_time.to_csv("SVM_reject_test_time.csv", index = False)


run_SVMrej(argv[1], argv[2], argv[3], argv[4], argv[5])
