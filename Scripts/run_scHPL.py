import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import time as tm
from scHPL import progressive_learning, predict, evaluate

def run_scHPL(DataPath, LabelsPath, TestPaths, OutputDirs):
	'''
	run scHPL
	Wrapper script to run scHPL on a benchmark dataset with 5-fold cross validation,
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv).
	OutputDir : Output directory defining the path of the predictions

	'''
	print("read the expression file")
	# read the data
	data = pd.read_csv(DataPath,index_col=0,sep=',')
	labels = pd.read_csv(LabelsPath, index_col = 0)
	# Map gene names to Upper 
	data.columns = map(str.upper, data.columns)
	# subset only to common genes
	# this removes dups
	data = data.groupby(level=0, axis=1).first()
	# Add assertion statement
	# get the Y vector (targets)
	y_train = np.array(labels['label'])
	train_time = []
	start = tm.time()
	classifier = 'svm'
	dimred = False
	threshold = 0.25
	tree = progressive_learning.learn_tree([data], [y_train], classifier = classifier, dimred = dimred, threshold = threshold)
	train_time.append(tm.time()-start)
	i = 0
	for test in TestPaths:
		test_time = []
		truelab = []
		pred = []
		OutputDir = OutputDirs[i]
		test = pd.read_csv(test, index_col = 0)
		test = test.groupby(level=0, axis=1).first()
		print(test.index)
		# predict
		start = tm.time()
		ypred = predict.predict_labels(test, tree)
		test_time.append(tm.time()-start)
		print(pred)
		prediction = pd.DataFrame({'scHPL':ypred})
		prediction.index = test.index
		tr_time = pd.DataFrame(train_time)
		ts_time = pd.DataFrame(test_time)
		path = os.path.join(OutputDir, "scHPL")
		os.chdir(path)
		prediction.to_csv("scHPL_pred.csv",index = True)
		tr_time.to_csv("scHPL_training_time.csv",index = False)
		ts_time.to_csv("scHPL_test_time.csv",index = False)
		i = i +1

  
parser = argparse.ArgumentParser()
parser.add_argument("-ref", "--ref", type=str, help="input reference path")
parser.add_argument("-labs", "--labs", type=str, help="input reference labels path")
parser.add_argument("-test", "--test", nargs="*", type=str, help="query expression matrices paths")
parser.add_argument("-output_dir", "--output_dir", nargs="*", type=str, help="output dirs for savig predictions")

args = parser.parse_args()
args.output_dir = [arg for arg in args.output_dir if os.path.isdir(arg) ]

# Function 
print(args.test)
print(args.output_dir)
run_scHPL(args.ref, args.labs, args.test, args.output_dir)
		
		
		






