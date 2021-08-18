# Single-Cell-Prediction-Pipeline


scRNA seq based Prediction of cell-types:

# How to use:

You need to have Snakemake installed in your linux enivornment. You need to also install all the dependancies for the tools you plan on using. You have to copy everything present in this repository and not break paths because it would disrupt the dependancies. One note is to change the paths in run_ACTINN.py to match your own directories when you clone the repository. The paths are in lines 44,45,49.

R Version 4.0.5 and Python 3.6.5 are used.

Current version of snakemake is snakemake/5.32.0

```bash
snakemake --use-conda --configfile config.yml --cores 3
```

# Config File:
```yaml 
output_dir: <path to outputs directory>
reference: <path to reference csv file with counts per cell, genes as columns and cells as rows>
labfile: <csv with labels per cell, the column header for the labels should be "label">
test: <path to test csv file with counts per cell, genes as columns and cells as rows>
rejection: <whether or not to reject poorly classified cells by SVM, default is True>
tools_to_run: # List of tools to run
  - <tool 1>
  - <tool 2>
  - <...>
```

## An Example Config is attached 

# Tools Available

1. ACTINN
2. SciBet
3. Spearman Correlation
4. SVM
5. SVM Rejection
6. SingleR
7. SingleCellNet

and many tools such as scMap Cell and my own classifier are being tested to be integrated in the pipeline.

# Submission File:

An example of the submission file is also available in this repository and is called submit.sh

# Packages Required:

## Python Specific Libraries:

```
tensorboard==1.7.0
tensorboard-data-server==0.6.1
tensorboard-plugin-wit==1.8.0
tensorflow==1.7.0
tensorflow-estimator==2.5.0
sklearn==0.0
scikit-learn==0.24.1
pandas==1.1.5
numpy==1.19.5
numpy-groupies==0.9.13
numpydoc==1.1.0
```

## R Libraries:

```
singleCellNet == 0.1.0
scibet == 1.0
SingleR == 1.4.1
Seurat == 4.0.3
dplyr == 1.0.7
tidyr == 1.1.3
viridis == 0.6.1
ggsci == 2.9
tidyverse == 1.3.1
```
