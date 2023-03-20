import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import time as tm
from scHPL import progressive_learning, predict, evaluate

def run_scHPL(RefPath, LabelsPath, QueryPaths, OutputDirs):
	"""
	Author: Hussein Lakkis
	Date: 2023-03-19
	run  classifier: scHPL 
	Wrapper script to run an scHPL classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
	Parameters
	----------
	RefPath : Reference Dataset file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell type annotations file path (.csv) with 'label' header .
	QueryPaths : Query dataset paths : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDirs : Output directory defining the path of the exported files for each query.
  	"""
	
	print("@Reading the reference expression file")

	# read the data
	Ref = pd.read_csv(RefPath,index_col=0,sep=',')
	labels = pd.read_csv(LabelsPath, index_col = 0)

	# Map gene names to Upper 
	Ref.columns = map(str.upper, Ref.columns)

	# this removes dups
	Ref = Ref.groupby(level=0, axis=1).first()

	# get the Y vector (targets)
	y_train = np.array(labels['label'])

	# set parameters for the classifier
	train_time = []
	start = tm.time()
	classifier = 'svm'
	dimred = False
	threshold = 0.25

	# train scHPL tree
	tree = progressive_learning.learn_tree([Ref], [y_train], classifier = classifier, dimred = dimred, threshold = threshold)
	train_time.append(tm.time()-start)

	# loop over query datasets
	i = 0

	for query in QueryPaths:
		query_time = []
		truelab = []
		pred = []
		
		# set output dir for each query sample
		OutputDir = OutputDirs[i]

		# read query
		query = pd.read_csv(query, index_col = 0)
		query = query.groupby(level=0, axis=1).first()
		print(query.index)

		# predict
		start = tm.time()
		ypred = predict.predict_labels(query, tree)
		query_time.append(tm.time()-start)

		# create prediction dataframe in addition to training and testing times
		prediction = pd.DataFrame({'scHPL' : ypred})
		prediction.index = query.index
		tr_time = pd.DataFrame(train_time)
		ts_time = pd.DataFrame(query_time)

		# write down output files
		path = os.path.join(OutputDir, "scHPL")
		os.chdir(path)
		prediction.to_csv("scHPL_pred.csv",index = True)
		tr_time.to_csv("scHPL_training_time.csv",index = False)
		ts_time.to_csv("scHPL_query_time.csv",index = False)
		i = i +1

# get command line args for this script
parser = argparse.ArgumentParser()
parser.add_argument("-ref", "--ref", type=str, help="input reference path")
parser.add_argument("-labs", "--labs", type=str, help="input reference labels path")
parser.add_argument("-query", "--query", nargs="*", type=str, help="query expression matrices paths")
parser.add_argument("-output_dir", "--output_dir", nargs="*", type=str, help="output dirs for savig predictions")

args = parser.parse_args()
args.output_dir = [arg for arg in args.output_dir if os.path.isdir(arg) ]

# execute the scHPL function
print(args.query)
print(args.output_dir)
run_scHPL(args.ref, args.labs, args.query, args.output_dir)
		
		
		






