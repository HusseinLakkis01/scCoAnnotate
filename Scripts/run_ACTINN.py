import os 
import numpy as np
import pandas as pd
import time as tm
from sys import argv
from pathlib import Path

def run_ACTINN(RefPath, LabelsPath, TestPath, OutputDir):
	'''
	run ACTINN
	Wrapper script to run ACTINN on a a reference and test datasets
	outputs list of predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	Reference : reference file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv).
	test : test file path (.csv), cells-genes matrix with cell unique barcodes 
	OutputDir : Output directory defining the path of the exported predictions

	'''
	# read the data sets
	reference = pd.read_csv(RefPath,index_col=0,sep=',')
	labels = pd.read_csv(LabelsPath, header = 0, index_col=0)
	test = pd.read_csv(TestPath,index_col=0)
	# folder with results
	os.chdir(OutputDir)
	tot=[]
	# transpose the datasets as ACTINN requires genes as rows
	train = reference.transpose()
	test = test.transpose()
	y_train = labels['label']
	# save files
	train.to_csv("train.csv")
	test.to_csv("test.csv")
	y_train.to_csv("train_lab.csv", header = False, index = True, sep = '\t')
	# remove these objects to save memory
	del(reference)
	del(train)
	del(test)
	tm.sleep(60)
	# RUN ACTINN formatting
	os.system("python /project/kleinman/hussein.lakkis/from_hydra/2021_07_15-Pipeline/Scripts/ACTINN_scripts/actinn_format.py -i train.csv -o train -f csv")   
	os.system("python /project/kleinman/hussein.lakkis/from_hydra/2021_07_15-Pipeline/Scripts/ACTINN_scripts/actinn_format.py -i test.csv -o test -f csv")
	# measure total time
	start = tm.time()
	# execute the actinn prediction file
	os.system("python /project/kleinman/hussein.lakkis/from_hydra/2021_07_15-Pipeline/Scripts/ACTINN_scripts/actinn_predict.py -trs train.h5 -trl train_lab.csv -ts test.h5")	
	tot.append(tm.time()-start)
	tm.sleep(60)
	# read predictions and probabilities
	predlabels = pd.read_csv('predicted_label.txt',header=0,index_col=0, sep='\t')
	prob = pd.read_csv('predicted_probabilities.txt',header=0,index_col=0, sep='\t')
	prob = prob.max()   
	pred = pd.DataFrame({"ACTINN_Prediction":predlabels['celltype'], "ACTINN_score":prob}, index = predlabels.index)
	tot_time = pd.DataFrame(tot)
	# output the files
	pred.to_csv("ACTINN_pred.csv", index = True)
	tot_time.to_csv("ACTINN_training_time.csv", index = False)
	tot_time.to_csv("ACTINN_test_time.csv", index = False)

# run function with arguments provided 
run_ACTINN(argv[1], argv[2], argv[3], argv[4])
# remove intermediate files to save space
os.remove("train.csv")
os.remove("test.csv")
os.remove("train.h5")
os.remove("train_lab.csv")
os.remove("test.h5")
os.remove("predicted_label.txt")
os.remove("predicted_probabilities.txt")