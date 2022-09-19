# import libraries
report: "report/workflow.rst"
import os

# Get the names of query samples from the paths given in the query section of the config
samples = [os.path.basename(os.path.dirname(query_path)) for query_path in config['query_datasets']]

# Get the names of tools to run from the list provided in the config
tools = config['tools_to_run']

# Create subdirectories for each query sample 
for sample in samples:
  path = "{output_dir}/{sample}".format(output_dir = config["output_dir"], sample = sample)
  if not os.path.exists(path):
    os.mkdir(path)

# create the prediction output list as an expand() output
output = expand(
  "{output_dir}/{sample}/{tool}/{tool}_pred.csv",
  tool = config['tools_to_run'],
  output_dir = config["output_dir"],
  sample  = samples)

# loop to only mark tools that havent been run on all samples to avoid rerunning
to_run = []
# loop over all tools
for tool in tools:
  alreadyrun = True
  for sample in samples:
    # check if the prediction output of the tool is present for all query datasets
    path = "{output_dir}/{sample}/{tool}/{tool}_pred.csv".format(output_dir = config["output_dir"], sample = sample, tool = tool)
    # if file is not present in all, mark it for running
    if not os.path.isfile(path):
      alreadyrun = False
  if alreadyrun == False:
    to_run.append(tool)

# config file to check if certain genes required by the user are present inb the list of common genes between reference and queries
genes_required = config['genes_required']

# create benchmarking directory is user chooses to perform benchmarking on the reference
benchmarking_dir = "{output_dir}/benchmarking".format(output_dir = config["output_dir"])
if not os.path.exists(benchmarking_dir) and bool(config['benchmark']):
  os.mkdir(benchmarking_dir)

# create the benchmarking config file and save to benchmarking subdirectory
if bool(config['benchmark']):
  outF = open("{dir}/benchmark.yml".format(dir = benchmarking_dir), "w")
  outF.write("output_dir: {output_dir}\n".format(output_dir = benchmarking_dir ))
  outF.write("datafile: {reference} \n".format(reference = config['training_reference']))
  outF.write("labfile: {labels}\n".format(labels = config['reference_annotations']))
  outF.write("rejection: \"True\"\n")
  outF.write("column: 2\n")
  outF.write("metrics: \n")
  outF.write("      - N1\n")
  outF.write("      - F2\n")
  outF.write("tools_to_run:\n")
  outF.write("      - Correlation\n")
  outF.write("      - SciBet\n")
  outF.close()

# invoke the benchmarking subworkflow (saved in the benchmarking folder in scCoAnnotate dir)
subworkflow Benchmark:
    workdir:
        "Benchmarking/"
    snakefile:
        "Benchmarking/Snakefile"
    configfile:
        "{benchmarking_dir}/benchmark.yml".format(benchmarking_dir =  benchmarking_dir)



# use rule all to dictate the final output
"""
One rule that directs the DAG which is represented in the rulegraph
"""
rule all:
  input:
      output = Benchmark("{benchmarking_dir}/evaluation/finished.csv".format(benchmarking_dir =  benchmarking_dir))
       if bool(config['benchmark']) else [],
      summary = expand(
        "{output_dir}/{sample}/Prediction_Summary.tsv",
        output_dir=config["output_dir"],
        sample  = samples),
      predictions = report(expand(
        "{output_dir}/{sample}/{tool}/{tool}_pred.csv",
        tool=to_run,
        output_dir=config["output_dir"],
        sample  = samples)),
      plot_and_embeddings = expand("{output_dir}/{sample}/figures/{tools}.png", 
        sample = samples, 
        output_dir = config['output_dir'], tools = config['tools_to_run']) if bool(config['plots']) else []


"""
 rule that gets the Consensus and concat all predictions
"""
rule plot:
  input:
      summary = expand(
        "{output_dir}/{sample}/Prediction_Summary.tsv",
        output_dir  = config["output_dir"],
        sample  = samples),      
      query =  expand("{query}", query=config['query_datasets']),
      output_dir =  expand("{output_dir}/{sample}", sample = samples, output_dir=config['output_dir'])
  output:
      expand("{output_dir}/{sample}/figures/{tools}.png", 
        sample = samples, 
        output_dir = config['output_dir'], tools = config['tools_to_run'])
    
  log: expand("{output_dir}/{sample}/plots.log", output_dir=config["output_dir"], sample  = samples)
  shell:
    "Rscript Scripts/plot_preds.R "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"    


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
  params:
    consensus = config.get("consensus")
  output:
    report(expand("{output_dir}/{sample}/Prediction_Summary.tsv", 
                    sample = samples, 
                    output_dir = config['output_dir']), caption = "report/summaries.rst", category = "Predictions")
  log: expand("{output_dir}/{sample}/Gatherpreds.log", output_dir=config["output_dir"], sample  = samples)
  shell:
    "python3 Scripts/Gather_Preds.py "
    "-i {input.sample} "
    "-c {params.consensus} "
    "-k {input.sample} "
    "&> {log}"    

"""
 rule that gets the gets the interesction in genes between samples and reference
 It outputs temporary reference and query datasets based on the common genes
"""
rule preprocess:
  input:
    reference = config['training_reference'],
    query =  expand("{query}", query=config['query_datasets']),
    output_dir = config['output_dir'],
  output:
    reference_result = temp("{output_dir}/expression.csv".format(output_dir = config['output_dir'])),
    query = temp(expand("{output_dir}/{sample}/expression.csv", output_dir=config["output_dir"],
        sample  = samples))
  params:
    check_genes = bool(config['check_genes']),
    genes_required = genes_required
  log: "{output_dir}/preprocess.log".format(output_dir=config['output_dir'])
  priority: 50
  shell:
    "Rscript Scripts/preprocess.R "
    "--ref {input.reference} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "--check_genes {params.check_genes} "
    "--genes_required {params.genes_required} "
    "&> {log}"


"""
Rules for R based tools.
"""
rule correlation:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/correlation/correlation_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/correlation/correlation_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/correlation/correlation_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/correlation/correlation.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_correlation.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"

rule SciBet:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/SciBet/SciBet_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/SciBet/SciBet_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SciBet/SciBet_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SciBet/SciBet.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SciBet.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"

rule scClassify:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scClassify/scClassify_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scClassify/scClassify_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scClassify/scClassify_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/scClassify/scClassify.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scClassify.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"
 
rule SingleCellNet:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  output:
    pred = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SingleCellNet/SingleCellNet.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SingleCellNet.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"   

rule scPred:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scPred/scPred_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scPred/scPred_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scPred/scPred_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/scPred/scPred.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scPred.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"   

rule SingleR:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/SingleR/SingleR_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/SingleR/SingleR_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SingleR/SingleR_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/SingleR/SingleR.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_SingleR.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"
    
rule CHETAH:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/CHETAH/CHETAH_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    total_time = expand("{output_dir}/{sample}/CHETAH/CHETAH_total_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/CHETAH/CHETAH.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_CHETAH.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"
    
rule scmapcluster:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scmapcluster/scmapcluster_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/scmapcluster/scmapcluster.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scmapcluster.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"

rule scmapcell:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/scmapcell/scmapcell_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scmapcell/scmapcell_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scmapcell/scmapcell_training_time.csv",sample  = samples,output_dir=config["output_dir"]) 
  
  log: expand("{output_dir}/{sample}/scmapcell/scmapcell.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "Rscript Scripts/run_scmapcell.R "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"

"""
Rules for python based tools.
"""

rule ACTINN:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])

  output:
    pred = expand("{output_dir}/{sample}/ACTINN/ACTINN_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/ACTINN/ACTINN_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/ACTINN/ACTINN_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/ACTINN/ACTINN.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_ACTINN.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}" 
    
rule SVM_reject:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  params:
    rejection = config.get("rejection", "")
  output:
    pred = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/SVM_reject/SVM_reject_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/SVM_reject/SVM_reject.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_SVM_reject.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
        "--rej {params.rejection}   "
    "--output_dir {input.output_dir} "
    "&> {log}"  

rule scHPL:
  input:
    reference = "{output_dir}/expression.csv".format(output_dir =config['output_dir']),
    labfile = config['reference_annotations'],
    query = expand("{output_dir}/{sample}/expression.csv",sample = samples,output_dir=config['output_dir']),
    output_dir =  expand("{output_dir}/{sample}",sample = samples,output_dir=config['output_dir'])
  output:
    pred = expand("{output_dir}/{sample}/scHPL/scHPL_pred.csv", sample  = samples,output_dir=config["output_dir"]),
    query_time = expand("{output_dir}/{sample}/scHPL/scHPL_query_time.csv",sample  = samples,output_dir=config["output_dir"]),
    training_time = expand("{output_dir}/{sample}/scHPL/scHPL_training_time.csv",sample  = samples,output_dir=config["output_dir"])
  log: expand("{output_dir}/{sample}/scHPL/scHPL.log", sample = samples,output_dir=config["output_dir"])
  shell:
    "python3 Scripts/run_scHPL.py "
    "--ref {input.reference} "
    "--labs {input.labfile} "
    "--query {input.query} "
    "--output_dir {input.output_dir} "
    "&> {log}"  
