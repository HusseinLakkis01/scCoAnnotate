# Single-Cell-Prediction-Pipeline <img src ="https://user-images.githubusercontent.com/59002771/130340419-3d1eff0b-ecb2-4104-9bf4-1bb968aff433.png" width="50" height="50">

scRNA seq based Prediction of cell-types using a fast and efficient pipeline to increase automation and reduce the need to run several scripts and experiments. The pipeline allows the user to select what single-cell projection tools they want to run.
The pipeline also features parallelization options to exploit the computational resources available. 

# Installation and Dependencies

Install Snakemake in your linux environment.

You need to have have R Version 4.0.5 and Python 3.6.5.

```bash
$ conda activate base
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```



You need to also install all the dependancies for the tools you plan on using. You have to copy everything present in this repository and not break paths because it would disrupt the dependancies. One note is to change the paths in run_ACTINN.py to match your own directories when you clone the repository. The paths are in lines 44,45,49.



Current version of snakemake is snakemake/5.32.0

# Quickstart

Using snakemake is straight forward and simple. The rules and processes are arranged as per this rule graph:

<img width="758" align = 'center' alt="rule_graph" src="https://user-images.githubusercontent.com/59002771/130340625-1239a7ec-dfd5-4005-aa90-c65ada201886.png">



You need to set everything up in a config file and then run the following command:

```bash
snakemake --use-conda --configfile config.yml --cores 3
```

##  Config File:
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

### An Example Config is attached 

```yaml 

output_dir: Results<img width="507" alt="Screen Shot 2021-08-21 at 10 52 30 PM" src="https://user-images.githubusercontent.com/59002771/130340335-ede960fd-05d2-4983-84a4-37bcd88743a1.png">
<img width="507" alt="Screen Shot 2021-08-21 at 10 52 30 PM" src="https://user-images.githubusercontent.com/59002771/130340339-ef8e027e-5b38-47a2-aaae-ff626e00e505.png">
![cell](https://user-images.githubusercontent.com/59002771/130340418-c4723b23-ea88-4cf2-ac73-97756bfe044f.png)

reference: /project/kleinman/hussein.lakkis/from_hydra/2021_01_07-Cross_Validation_and_Benchmark/2021_04_05-SVM_and_SVMrej/data/scRNAseq_Benchmark_datasets/Joint_Mouse/joint_mouse.training.csv
labfile: /project/kleinman/hussein.lakkis/from_hydra/2021_01_07-Cross_Validation_and_Benchmark/2021_04_05-SVM_and_SVMrej/data/scRNAseq_Benchmark_datasets/Joint_Mouse/full_labels.csv
test: /project/kleinman/zahedeh.bashardanesh/from_beluga/2020-11_MANAV/data/S-10068_28741/expr.csv
rejection: "True"
tools_to_run:
      - correlation
      - SciBet
      - ACTINN
      - SVM_reject
```

## Submission File:

An example of the submission file is also available in this repository and is called submit.sh

``` bash 
#!/usr/bin/bash
#PBS -N Snakemake)_Pipeline
#PBS -o logs/err.txt
#PBS -e logs/out.txt
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=3
#PBS -l mem=125G
#PBS -l vmem=125G

# you need to do this step 

cd {root directory of the pipeline}
mkdir -p logs

# set up the environment
#conda init bash
module load miniconda/3.8
source ~/.conda_init.sh
module load snakemake/5.32.0
module load hydra/R/4.0.5
module load python/3.6.5

# Run the snakemake
snakemake --use-conda --configfile config.yml --cores 3
```

# Tools Available

1. ACTINN
2. SciBet
3. Spearman Correlation
4. SVM
5. SVM Rejection
6. SingleR
7. SingleCellNet

and many tools such as scMap Cell and my own classifier are being tested to be integrated in the pipeline.



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
# Adding New Tools:
