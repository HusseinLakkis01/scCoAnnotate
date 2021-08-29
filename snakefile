"""
One rule that directs the DAG which is represented in the rulegraph
"""
rule all:
  input:
      expand(
        "{output_dir}/Prediction_Summary.tsv",
        output_dir=config["output_dir"]
)
"""
 rule that gets the Consensus 
"""
rule concat:
  input:
      expand(
        "{output_dir}/{tool}/{tool}_pred.csv",
        tool=config["tools_to_run"],
        output_dir=config["output_dir"])
  output:
    result = "{output_dir}/Prediction_Summary.tsv"
  log: "{output_dir}/Gatherpreds.log"
  shell:
    "python3 Gather_Preds.py "
    "{wildcards.output_dir}"
    "&> {log}"    

"""
Rule for R based tools.
"""
rule Correlation:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"]
  output:
    pred = "{output_dir}/correlation/correlation_pred.csv",
    test_time = "{output_dir}/correlation/correlation_test_time.csv",
    training_time = "{output_dir}/correlation/correlation_training_time.csv"
  log: "{output_dir}/correlation/correlation.log"
  shell:
    "Rscript Scripts/run_Correlation.R "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/correlation "
    "&> {log}"

rule SciBet:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"]
  output:
    pred = "{output_dir}/SciBet/SciBet_pred.csv",
    test_time = "{output_dir}/SciBet/SciBet_test_time.csv",
    training_time = "{output_dir}/SciBet/SciBet_training_time.csv"
  log: "{output_dir}/SciBet/SciBet.log"
  shell:
    "Rscript Scripts/run_SciBet.R "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/SciBet "
    "&> {log}"

rule ACTINN:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"]
  output:
    pred = "{output_dir}/ACTINN/ACTINN_pred.csv",
    test_time = "{output_dir}/ACTINN/ACTINN_test_time.csv",
    training_time = "{output_dir}/ACTINN/ACTINN_training_time.csv"
  log: "{output_dir}/ACTINN/ACTINN.log"
  shell:
    "python3 Scripts/run_ACTINN.py "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/ACTINN "
    "&> {log}"

rule SVM:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"],
  output:
    pred = "{output_dir}/SVM/SVM_pred.csv",
    test_time = "{output_dir}/SVM/SVM_test_time.csv",
    training_time = "{output_dir}/SVM/SVM_training_time.csv"
  log: "{output_dir}/SVM/SVM.log"
  shell:
    "python3 Scripts/run_SVM.py "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/SVM "
    "&> {log}"

rule SVM_reject:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"],
  params:
    rejection = config.get("rejection", "")
  output:
    pred = "{output_dir}/SVM_reject/SVM_reject_pred.csv",
    test_time = "{output_dir}/SVM_reject/SVM_reject_test_time.csv",
    training_time = "{output_dir}/SVM_reject/SVM_reject_training_time.csv"
  log: "{output_dir}/SVM_reject/SVM_reject.log"
  shell:
    "python3 Scripts/run_SVM_reject.py "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/SVM_reject "
    "{params.rejection} "
    "&> {log}"

rule singleCellNet:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"]
  output:
    pred = "{output_dir}/singleCellNet/singleCellNet_pred.csv",
    test_time = "{output_dir}/singleCellNet/singleCellNet_test_time.csv",
    training_time = "{output_dir}/singleCellNet/singleCellNet_training_time.csv"
  log: "{output_dir}/singleCellNet/singleCellNet.log"
  shell:
    "Rscript Scripts/run_singleCellNet.R "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/singleCellNet "
    "&> {log}"

rule scPred:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"]
  output:
    pred = "{output_dir}/scPred/scPred_pred.csv",
    test_time = "{output_dir}/scPred/scPred_test_time.csv",
    training_time = "{output_dir}/scPred/scPred_training_time.csv"
  log: "{output_dir}/scPred/scPred.log"
  shell:
    "Rscript Scripts/run_scPred.R "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/scPred "
    "&> {log}"

rule SingleR:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"]
  output:
    pred = "{output_dir}/SingleR/SingleR_pred.csv",
    test_time = "{output_dir}/SingleR/SingleR_test_time.csv",
    training_time = "{output_dir}/SingleR/SingleR_training_time.csv"
  log: "{output_dir}/SingleR/SingleR.log"
  shell:
    "Rscript Scripts/run_SingleR.R "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/SingleR "
    "&> {log}"
    
    
rule SingleCellNet:
  input:
    reference = config["reference"],
    labfile = config["labfile"],
    test =  config["test"]
  output:
    pred = "{output_dir}/SingleCellNet/SingleCellNet_pred.csv",
    test_time = "{output_dir}/SingleCellNet/SingleCellNet_test_time.csv",
    training_time = "{output_dir}/SingleCellNet/SingleCellNet_training_time.csv"
  log: "{output_dir}/SingleCellNet/SingleCellNet.log"
  shell:
    "Rscript Scripts/run_SingleCellNet.R "
    "{input.reference} "
    "{input.labfile} "
    "{input.test} "
    "{wildcards.output_dir}/SingleCellNet "
    "&> {log}"
