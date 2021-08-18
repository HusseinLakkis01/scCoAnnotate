import os 
import numpy as np
import pandas as pd
from sys import argv
from pathlib import Path
import re

# Function 
def glob_re(path, regex="_pred.csv", glob_mask="**/*", inverse=False):
	'''
	This is a function that looks for all files in the given directory and sub dirs that
	match a specific regex pattern
	The goal is to get prediction files
	'''
	p = Path(path)
	if inverse:
		res = [str(f) for f in p.glob(glob_mask) if not re.search(regex, str(f))]
	else:
		res = [str(f) for f in p.glob(glob_mask) if re.search(regex, str(f))]
	return res



def run_concat(Results_Folder_Path):
	'''
	run Concat
	Wrapper to collect all predictions
	outputs lists of true and predicted cell labels as csv files, with all tools
  
	Parameters
	----------
	Results_Folder_Path: Path to results folder
	'''
	paths = glob_re(Results_Folder_Path, regex="_pred.csv", glob_mask="**/*", inverse=False)
	result = pd.read_csv(paths[0],index_col = 0)
	for path in paths[1:]:
		current = pd.read_csv(path, index_col = 0)
		result = pd.concat([result, current], axis = 1)
	mode = result.select_dtypes(exclude=['float64']).mode(axis = 1)
	result['Consensus'] = np.where(mode.isna().any(1), mode[0], 'Unsure')
	os.chdir(Results_Folder_Path)
	result.to_csv("Prediction_Summary.tsv", sep = "\t", index = True)


run_concat(argv[1])
