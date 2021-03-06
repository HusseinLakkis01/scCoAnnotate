import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
import time as tm

def run_SVM_reject(RefPath, LabelsPath, TestPaths, OutputDirs, rejection):
	'''
	run SVM_reject
	Wrapper script to run SVM_reject on a benchmark dataset with 5-fold cross validation,
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv).
	OutputDir : Output directory defining the path of the predictions

	'''
	# get rejectiom
	reject = bool(rejection)
	threshold = 0.6
	# read the data
	print("read the expression file")
	# read the data
	data = pd.read_csv(RefPath,index_col=0,sep=',')
	labels = pd.read_csv(LabelsPath, index_col = 0)
	# Map gene names to Upper 
	data.columns = map(str.upper, data.columns)
	# subset only to common genes
	# this removes dups
	data = data.groupby(level=0, axis=1).first()
	data = np.log1p(data)
	# get the Y vector (targets)
	y_train = np.array(labels['label'])
	train_time=[]
	clf = LinearSVC()
	Classifier = CalibratedClassifierCV(clf)
	start=tm.time()
	# fit and train
	Classifier.fit(data, y_train)
	# get training time
	train_time.append(tm.time()-start)
	start=tm.time()

	i = 0
	for test in TestPaths:
		test_time = []
		OutputDir = OutputDirs[i]
		test = pd.read_csv(test, index_col = 0)
		test = test.groupby(level=0, axis=1).first()
		test = np.log1p(test)
		# predict
		# get preds
		predicted = Classifier.predict(test)
		# predict the probabilities 
		prob = Classifier.predict_proba(test)
		test_time.append(tm.time()-start)
		probabilities = np.max(prob, axis = 1)
		# create a dataframe with the output
		if reject:
			# remove bad cells
			unlabeled = np.where(probabilities < threshold)
			predicted[unlabeled] = 'Unknown'
			pred = pd.DataFrame({'SVM_Predictions_Rej' : predicted,'SVM_Predictions':Classifier.predict(test), 'SVM_Probability_score' : probabilities})
			pred.index = test.index
		else:
			pred = pd.DataFrame({'SVM_Predictions' : predicted, 'SVM_Probability_score' : probabilities})
		
		test_time.append(tm.time()-start)
		pred.index = test.index
		train_time = pd.DataFrame(train_time)
		ts_time = pd.DataFrame(test_time)
		path = os.path.join(OutputDir, "SVM_reject")
		os.chdir(path)
		pred.to_csv("SVM_reject_pred.csv",index = True)
		train_time.to_csv("SVM_reject_training_time.csv",index = False)
		ts_time.to_csv("SVM_reject_test_time.csv",index = False)
		i = i +1

  
parser = argparse.ArgumentParser()
parser.add_argument("-ref", "--ref", type=str, help="input reference path")
parser.add_argument("-labs", "--labs", type=str, help="input reference labels path")
parser.add_argument("-test", "--test", nargs="*", type=str, help="query expression matrices paths")
parser.add_argument("-output_dir", "--output_dir", nargs="*", type=str, help="output dirs for savig predictions")
parser.add_argument("-rej", "--rejection", type=bool, help="output dirs for savig predictions")

args = parser.parse_args()
args.output_dir = [arg for arg in args.output_dir if os.path.isdir(arg) ]

# Function 
print(args.test)
print(args.output_dir)
run_SVM_reject(args.ref, args.labs, args.test, args.output_dir, args.rejection)






