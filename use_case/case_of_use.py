from lib.dorotheadb import *
import random
import os
import pandas
import matplotlib.pyplot as plt
import glob
import numpy
import json
from collections import defaultdict
import multiprocessing
from time import process_time

# Select the annotation to use
annotations = ["KEGG","GO_BP","Reactome","WikiPathways"]

for annotation in annotations:
  print(annotation)
  outFile = "case_of_use/SLE_"+annotation+".tsv"
  os.system("python3 TFsEnrichment.py --tfs tfs_SLE.txt --output "+outFile+" --annotation "+annotation)
  outFile = "case_of_use/Cancer_"+annotation+".tsv"
  os.system("python3 TFsEnrichment.py --tfs tfs_Cancer.txt --output "+outFile+" --annotation "+annotation)
