import os
# Get the names of query samples from the paths given in the test section of the config
samples = [os.path.basename(os.path.dirname(test_path)) for test_path in config["test"]]
tools = config['tools_to_run']
# Create subdirectories for each query sample
for sample in samples:
  path = "{output_dir}/{sample}".format(output_dir = config["output_dir"], sample = sample)
  if not os.path.exists(path):
    os.mkdir(path)

output = expand(
  "{output_dir}/{sample}/{tool}/{tool}_pred.csv",
  tool= config['tools_to_run'],
  output_dir=config["output_dir"],
  sample  = samples)

to_run = []
for tool in tools:
  alreadyrun = True
  for sample in samples:
    path = "{output_dir}/{sample}/{tool}/{tool}_pred.csv".format(output_dir = config["output_dir"], sample = sample, tool = tool)
    if not os.path.isfile(path):
      alreadyrun = False
  if alreadyrun == False:
    to_run.append(tool)





"""
One rule that directs the DAG which is represented in the rulegraph
"""
rule all:
  input:
      res1 = expand(
        "{output_dir}/{sample}/Prediction_Summary.tsv",
        output_dir=config["output_dir"],
        sample  = samples),
      res2= expand(
        "{output_dir}/{sample}/{tool}/{tool}_pred.csv",
        tool=to_run,
        output_dir=config["output_dir"],
        sample  = samples)

"""
 rule that gets the Consensus 
"""
rule concat:
  input:
      results = expand("{output_dir}/{sample}/{tool}/{tool}_pred.csv",
        tool=to_run,
        output_dir=config["output_dir"],
        sample  = samples),
      sample =  expand("{output_dir}/{sample}",
        output_dir=config["output_dir"],
        sample  = samples)
  output:
    result = expand("{output_dir}/{sample}/Prediction_Summary.tsv", 
                    sample = samples, 
                    output_dir = config['output_dir'])
  log: expand("{output_dir}/{sample}/Gatherpreds.log", output_dir=config["output_dir"], sample  = samples)
  shell:
    "python3 Scripts/Gather_Preds.py "
    "-i {input.sample} "
    "&> {log}"    

"""
 rule that gets the gets the interesction in genes between samples and reference
 It outputs temporary reference and query datasets based on the common genes
"""
rule preprocess:
  input:
    reference = config["reference"],
    test =  expand("{test}", test=config["test"]),
    output_dir = config['output_dir']
  output:
    reference_result = temp("{output_dir}/expression.csv".format(output_dir = config['output_dir'])),
    test = temp(expand("{output_dir}/{sample}/expression.csv", output_dir=config["output_dir"],
        sample  = samples))
  log: "{output_dir}/preprocess.log".format(output_dir=config['output_dir'])
  priority: 50
  shell:
    "Rscript Scripts/preprocess.R "
    "--ref {input.reference} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"


"""
Rules for R based tools.
"""
rule Correlation:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/correlation/correlation_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/correlation/correlation_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/correlation/correlation_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/correlation/correlation.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_Correlation.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"

rule SciBet:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/SciBet/SciBet_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/SciBet/SciBet_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SciBet/SciBet_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SciBet/SciBet.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SciBet.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"

rule scClassify:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scClassify/scClassify_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/scClassify/scClassify_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scClassify/scClassify_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/scClassify/scClassify.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scClassify.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"
 
rule SingleCellNet:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  output:
    pred = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SingleCellNet.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"   

rule scPred:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scPred/scPred_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/scPred/scPred_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scPred/scPred_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/scPred/scPred.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scPred.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"   

rule SingleR:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/SingleR/SingleR_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/SingleR/SingleR_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SingleR/SingleR_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/SingleR/SingleR.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SingleR.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"
    
rule CHETAH:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/CHETAH/CHETAH_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    total_time = expand("{output_dir}/{sample}/CHETAH/CHETAH_total_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/CHETAH/CHETAH.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_CHETAH.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"
    
rule scmapcluster:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/scmapcluster/scmapcluster.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scmapcluster.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"

rule scmapcell:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scmapcell/scmapcell_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/scmapcell/scmapcell_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scmapcell/scmapcell_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/scmapcell/scmapcell.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scmapcell.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"

"""
Rules for python based tools.
"""

rule ACTINN:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/ACTINN/ACTINN_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/ACTINN/ACTINN_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/ACTINN/ACTINN_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/ACTINN/ACTINN.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_ACTINN.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}" 
    
rule SVM_reject:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  params:
    rejection = config.get("rejection", "")
  output:
    pred = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SVM_reject/SVM_reject.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_SVM_reject.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
        "--rej {params.rejection}   "
    "--output_dir {input.output_dir} "
    "&> {log}"  

rule scHPL:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config["labfile"],
    test = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  output:
    pred = expand("{output_dir}/{sample}/scHPL/scHPL_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    test_time = expand("{output_dir}/{sample}/scHPL/scHPL_test_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scHPL/scHPL_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/scHPL/scHPL.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_scHPL.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--test {input.test} "
    "--output_dir {input.output_dir} "
    "&> {log}"  