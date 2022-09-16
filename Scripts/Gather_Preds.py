import os 
import numpy as np
import pandas as pd
from pathlib import Path
import re
import argparse

# get arguments provided on the command line
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", nargs="*", type=str, help="input directory path")
parser.add_argument("-c", "--consensus", nargs="*", type=str, help="tools used for consensus", default = "all")
parser.add_argument("-k", "--splitter", nargs="*", type=str, help="string to account for variable input length")

# parse the arguments
args = parser.parse_args()

# Functions in this script

def glob_re(path, regex="_pred.csv", glob_mask="**/*", inverse=False):
	'''
	This is a function that looks for all files in the given directory and sub dirs that
	match a specific regex pattern
	The goal is to get prediction files
	'''
	p = Path(path)

	# get prediction files with regex (_pred.csv)
	if inverse:
		res = [str(f) for f in p.glob(glob_mask) if not re.search(regex, str(f))]
	else:
		res = [str(f) for f in p.glob(glob_mask) if re.search(regex, str(f))]
	return res



def run_concat(Results_Folder_Path, tools_for_consensus):
	'''
	run Concat
	Wrapper function to collect all predictions
	outputs lists of true and predicted cell labels as csv files, with all tools in addition to a consensus defined by user (majority vote)

	Parameters
	----------
	Results_Folder_Path: Path to results folder defined in the config
	tools_for_consensus: parameter defined by the user in the config to get a consensus prediction, default is to use all classifiers

	Output
	-------
	writes a tsv file containing all predictions from individuals tools in addition to the consensus predictions for all query samples

	'''
	# get prediction files per tools for a sample

	paths = glob_re(Results_Folder_Path, regex="_pred.csv", glob_mask="**/*", inverse=False)
	print(paths)
	print(paths[0])

	# append all predictions to one dataframe

	result = pd.read_csv(paths[0],index_col = 0)
	for path in paths[1:]:
		current = pd.read_csv(path, index_col = 0)
		result = pd.concat([result, current], axis = 1)

	# case 1: get consensus of all tools 
	if tools_for_consensus[0] == 'all':
		mode = result.select_dtypes(exclude=['float64']).mode(axis = 1)
		result['Consensus'] = np.where(mode.isna().any(1), mode[0], 'Unsure')

	# case 2: user defines a set of tools for consensus ( more than one)

	elif len(tools_for_consensus) > 1:

		# get consensus between these tools only
		consensus_df = pd.read_csv("{i}/{j[0]}/{j[0]}_pred.csv".format(i = Results_Folder_Path, j = tools_for_consensus), index_col = 0)

		for tool in tools_for_consensus[1:]:
			current = pd.read_csv("{i}/{j}/{j}_pred.csv".format(i = Results_Folder_Path, j =tool), index_col = 0)
			consensus_df = pd.concat([consensus_df, current], axis = 1)

			# find the majority vote
			mode = consensus_df.select_dtypes(exclude=['float64']).mode(axis = 1)
			print(mode)
			# append it to the results dataframe
			result['Consensus'] = np.where(mode.isna().any(1), mode[0], 'Unsure')
	
	# case 3: user accidenly enters one tool for consensus 
	else:
		result['Consensus'] = result.iloc[:, 0]

	# write results to each individual query subdir in the output directory
	os.chdir(Results_Folder_Path)
	result.index.name = "cellname"
	result.to_csv("Prediction_Summary.tsv", sep = "\t", index = True)


# run concat on the different samples
args.input = [arg for arg in args.input if os.path.isdir(arg) ]
for arg in args.input:
	print(arg)
	run_concat(arg,args.consensus)


