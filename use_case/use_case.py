import os

# Select the annotation to use
annotations = ["KEGG","GO_BP","Reactome","WikiPathways"]
methods = ["targetsF","targetsW"]
for method in methods:
  for annotation in annotations:
    print(annotation)
    outFile = "use_case/results/Cancer_"+annotation+"_"+method+".tsv"
    os.system("python3 TFsEnrichment.py -tfs use_case/tfs_Cancer.txt -m"+method+" -out "+outFile+" -ann "+annotation)
    outFile = "use_case/results/SLE_"+annotation+"_"+method+".tsv"
    os.system("python3 TFsEnrichment.py -tfs use_case/tfs_SLE.txt -m"+method+" -out "+outFile+" -ann "+annotation)
