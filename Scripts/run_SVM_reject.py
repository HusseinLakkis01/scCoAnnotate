import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
import time as tm

def run_SVM_reject(RefPath, LabelsPath, QueryPaths, OutputDirs, rejection):
	'''
	run SVM_reject
	Wrapper script to train SVM_reject on a referenxe dataset with 5-fold cross validation,
	outputs lists of predicted query cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath : Reference dataset file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv).
	QueryPaths: Query datasets file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDir : Output directory defining the path of the predictions.
	rejection: Option to reject cells with prediction probability below thresold: Bool

	'''
	# get rejectiom
	reject = bool(rejection)
	threshold = 0.6
	# read the reference
	print("read the expression file")
	# read the data
	ref = pd.read_csv(RefPath,index_col=0,sep=',')
	labels = pd.read_csv(LabelsPath, index_col = 0)
	# Map gene names to Upper 
	ref.columns = map(str.upper, ref.columns)
	# subset only to common genes
	# this removes dups
	ref = ref.groupby(level=0, axis=1).first()
	ref = np.log1p(ref)
	# get the Y vector (targets)
	y_train = np.array(labels['label'])
	train_time=[]
	clf = LinearSVC()
	Classifier = CalibratedClassifierCV(clf)
	start=tm.time()
	# fit and train
	Classifier.fit(ref, y_train)
	# get training time
	train_time.append(tm.time()-start)
	start=tm.time()

	i = 0
	for query in QueryPaths:
		query_time = []
		OutputDir = OutputDirs[i]
		query = pd.read_csv(query, index_col = 0)
		query = query.groupby(level=0, axis=1).first()
		query = np.log1p(query)
		# predict
		# get preds
		predicted = Classifier.predict(query)
		# predict the probabilities 
		prob = Classifier.predict_proba(query)
		query_time.append(tm.time()-start)
		probabilities = np.max(prob, axis = 1)
		# create a dataframe with the output
		if reject:
			# remove bad cells
			unlabeled = np.where(probabilities < threshold)
			predicted[unlabeled] = 'Unknown'
			pred = pd.DataFrame({'SVM_reject' : predicted,'SVM':Classifier.predict(query), 'SVM_Probability_score' : probabilities})
			pred.index = query.index
		else:
			pred = pd.DataFrame({'SVM' : predicted, 'SVM_Probability_score' : probabilities})
		
		query_time.append(tm.time()-start)
		pred.index = query.index
		train_time = pd.DataFrame(train_time)
		ts_time = pd.DataFrame(query_time)
		path = os.path.join(OutputDir, "SVM_reject")
		os.chdir(path)
		pred.to_csv("SVM_reject_pred.csv",index = True)
		train_time.to_csv("SVM_reject_training_time.csv",index = False)
		ts_time.to_csv("SVM_reject_query_time.csv",index = False)
		i = i +1

  
parser = argparse.ArgumentParser()
parser.add_argument("-ref", "--ref", type=str, help="input reference path")
parser.add_argument("-labs", "--labs", type=str, help="input reference labels path")
parser.add_argument("-query", "--query", nargs="*", type=str, help="query expression matrices paths")
parser.add_argument("-output_dir", "--output_dir", nargs="*", type=str, help="output dirs for savig predictions")
parser.add_argument("-rej", "--rejection", type=bool, help="output dirs for savig predictions", default = False)

args = parser.parse_args()
args.output_dir = [arg for arg in args.output_dir if os.path.isdir(arg) ]

# Function 
print(args.query)
print(args.output_dir)
run_SVM_reject(args.ref, args.labs, args.query, args.output_dir, args.rejection)






