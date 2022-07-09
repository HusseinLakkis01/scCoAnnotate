# scCoAnnotate <img src ="https://user-images.githubusercontent.com/59002771/130340419-3d1eff0b-ecb2-4104-9bf4-1bb968aff433.png" width="50" height="50">

# Summary

scRNA seq based Prediction of cell-types using a fast and efficient pipeline to increase automation and reduce the need to run several scripts and experiments. The pipeline allows the user to select what single-cell projection tools they want to run on a selected reference to annotate a list of query datasets. It then outputs a consensus of the predictions across tools selected. This pipeline trains classifiers on genes common to the reference and all query datasets. 

The pipeline also features parallelization options to exploit the computational resources available. 

# Installation and Dependencies

Install [Snakemake](https://snakemake.readthedocs.io/en/stable/) in your linux environment.

You need to have have [R](https://www.r-project.org/) Version 4.0.5 and Python 3.6.5.

```bash
$ conda activate base
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```



You need to also install all the dependancies for the tools you plan on using. You have to copy everything present in this repository and not break paths because it would disrupt the dependancies. One note is to change the paths in run_ACTINN.py to match your own directories when you clone the repository. The paths are in lines 44,45,49.



Current version of snakemake is snakemake/5.32.0

# Quickstart

Using snakemake is straight forward and simple. The rules and processes are arranged as per this rule graph:

Rule preprocess gets the common genes and creates temporary reference and query datasets based ob the common genes. Rule concat appends all predictions into one tab seperate file (prediction_summary.tsv) and gets the consensus prediction


![save_as_a_png](https://user-images.githubusercontent.com/59002771/178054140-e7129733-6a8f-4819-8162-c29b3954d303.png)




You need to set everything up in a config file and then run the following command:

```bash
snakemake --use-conda --configfile config.yml --cores 3
```

##  Config File:
```yaml 
output_dir: <path to outputs directory>
reference: <path to reference csv file with RAW counts per cell, genes as columns and cells as rows>
labfile: <csv with labels per cell, the column header for the labels should be "label">
test: - <path to test csv file 1 with RAW counts per cell, genes as columns and cells as rows>
      - <path to test csv file 2 with RAW counts per cell, genes as columns and cells as rows>
      .
      .
      .
rejection: <whether or not to reject poorly classified cells by SVM, default is True>
tools_to_run: # List of tools to run
  - <tool 1>
  - <tool 2>
  - <...>
```

### An Example Config is attached 

```yaml 
output_dir: /project/kleinman/hussein.lakkis/from_hydra/test
reference: /projects/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/reference/reference.csv
labfile: /projects/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/reference/labels.csv
test:
      - /projects/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/BT2016062/expression.csv
      - /projects/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/BT2018022/expression.csv
      - /projects/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/P-1190_S-1197/expression.csv
      - /projects/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/P-1569_S-1569/expression.csv
      - /projects/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/P-1694_S-1694_multiome/expression.csv
      - /projects/kleinman/hussein.lakkis/from_hydra/Collab/HGG_Selin_Revision/test/P-1701_S-1701_multiome/expression.csv
rejection: True
tools_to_run:
      - scmapcell
      - scmapcluster
      - ACTINN
      - SVM_reject
      - SingleCellNet
      - SciBet
      - scHPL
      - correlation
      - CHETAH
      - correlation
```

## Submission File:

An example of the submission file is also available in this repository and is called submit.sh

``` bash 
#!/usr/bin/bash
#PBS -N scCoAnnotate
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

1. [ACTINN](https://github.com/mafeiyang/ACTINN)
2. [SciBet](https://github.com/PaulingLiu/scibet)
4. [Spearman Correlation](https://statistics.laerd.com/statistical-guides/spearmans-rank-order-correlation-statistical-guide.php)
5. [SVM](https://scikit-learn.org/stable/modules/svm.html)
6. SVM Rejection
7. [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html)
8. [SingleCellNet](https://github.com/pcahan1/singleCellNet)
9. [CHETAH](https://www.bioconductor.org/packages/release/bioc/html/CHETAH.html)
10. [scHPL](https://github.com/lcmmichielsen/scHPL)
11. [scPred](https://github.com/powellgenomicslab/scPred)
12. [scmap (cell and cluster)](https://bioconductor.org/packages/release/bioc/html/scmap.html)

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
scHPL==0.0.2
```

## R Libraries:

```
scPred_1.9.2
SingleCellExperiment_1.12.0
SummarizedExperiment_1.20.0
CHETAH_1.6.0
scmap_1.12.0 
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

to add new tools, you have to add this template to the the snakefile as such:
``` python
rule {tool_name}:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/{tool_name}/{tool_name}_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/{tool_name}/{tool_name}_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/{tool_name}/{tool_name}_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/{tool_name}/{tool_name}.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_{tool_name}.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"
 ```   
 The tool script you add must generate outputs that match the output of the rule..



