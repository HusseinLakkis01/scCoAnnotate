import os 
import numpy as np
import pandas as pd
import time as tm
from sys import argv
from pathlib import Path
import argparse
import shutil


def run_ACTINN(RefPath, LabelsPath, QueryPaths, OutputDirs):
	"""
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: ACTINN 
	Wrapper script to run an ACTINN classifier 
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
	# create temp folder for ACTINN output
	tmp = "temp"
	tmp_path = os.path.abspath(tmp)
	if not os.path.exists(tmp_path):
		os.mkdir(tmp)

	# read the data sets
	reference = pd.read_csv(RefPath,index_col=0,sep=',')
	labels = pd.read_csv(LabelsPath, header = 0, index_col=0)
	os.chdir(tmp_path)

	# folder with results
	tot=[]
	reference = np.log1p(reference)
	# transpose the datasets as ACTINN requires genes as rows
	train = reference.transpose()
	y_train = labels['label']

	# save files to run ACTINN scripts
	train.to_csv("train.csv")
	y_train.to_csv("train_lab.csv", header = False, index = True, sep = '\t')

	# remove these objects to save memory
	del(reference)
	del(train)
	tm.sleep(60)

	# RUN ACTINN formatting
	os.system("python ../Scripts/ACTINN_scripts/actinn_format.py -i train.csv -o train -f csv")

	i = 0
	for Query in QueryPaths:
		print(Query) 
		OutputDir = OutputDirs[i]
		path = os.path.join(OutputDir, "ACTINN")

		# read and normalize
		os.chdir(tmp_path)
		Query = pd.read_csv(Query, index_col = 0)
		Query = np.log1p(Query)
		Query = Query.transpose()
		Query.to_csv("Query.csv")
		os.system("python ../Scripts/ACTINN_scripts/actinn_format.py -i Query.csv -o Query -f csv")
		# measure total time
		start = tm.time()

		# execute the actinn prediction file
		os.system("python ../Scripts/ACTINN_scripts/actinn_predict.py -trs train.h5 -trl train_lab.csv -ts Query.h5 -op False ")	
		tot.append(tm.time()-start)
		tm.sleep(60)
		
		# read predictions and probabilities
		predlabels = pd.read_csv('{tmp}/predicted_label.txt'.format(tmp = tmp_path),header=0,index_col=0, sep='\t')

		pred = pd.DataFrame({"ACTINN":predlabels['celltype']}, index = predlabels.index)
		tot_time = pd.DataFrame(tot)

		# output the files
		os.chdir(path)
		pred.to_csv("ACTINN_pred.csv", index = True)
		tot_time.to_csv("ACTINN_training_time.csv", index = False)
		tot_time.to_csv("ACTINN_query_time.csv", index = False)
		i = i+1


# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-ref", "--ref", type=str, help="input reference path")
parser.add_argument("-labs", "--labs", type=str, help="input reference labels path")
parser.add_argument("-query", "--query", nargs="*", type=str, help="query expression matrices paths")
parser.add_argument("-output_dir", "--output_dir", nargs="*", type=str, help="output dirs for savig predictions")

# get the parses
args = parser.parse_args()

# keep output dir with a full path
args.output_dir = [arg for arg in args.output_dir if os.path.isdir(arg) ]


# run function with arguments provided 
run_ACTINN(args.ref, args.labs, args.query, args.output_dir)

