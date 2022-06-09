import pandas
import os
from collections import Counter
from collections import defaultdict
import sys
import argparse
import numpy
from lib.dorotheadb import *
from lib.calculate_stats import *
from statsmodels.stats.multitest import multipletests
import json

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-db", "--database", help="select level of confidence for dorothea", default="C")
parser.add_argument("-org", "--organism", type=int, help="select tax id: 9606 for human and 10090 for mouse", default=9606)
parser.add_argument("-ann", "--annotation", type=str, help="select annotation to use", default="KEGG")
parser.add_argument("-out", "--output", type=str, help="output file", default="out.txt")
parser.add_argument("-tfs", "--tfs", type=str, help="list of tfs divided by commas (STAT1,STAT2,IRF1) or file to read in txt format")
args = parser.parse_args()
config = vars(args)

tfs = config["tfs"]
out = config["output"]
annotation = config["annotation"]
database = config["database"]
organism = config["organism"]

if tfs.endswith(".txt"):
    if os.path.exists(tfs):
        tfs = pandas.read_csv(tfs,header=None)[0].tolist()
        if type(tfs) != list:
            print("Error")
            print("Incorrect format file. One column with TFs names is required")
            sys.exit()
elif type(tfs.split(",")) == list:
    tfs = tfs.split(",")
else:
    print("Error")
    print("Incorrect format file. One column with TFs names is required")



# Read databases
dorotheadb = getDorotheaDB(confidence=database, org=organism)
annFileTable = "data/"+annotation+".tsv"
annTable = pandas.read_csv(annFileTable, sep="\t")
annTable = annTable.loc[annTable.organism == organism]

# Filter databases
dorotheadb = dorotheadb.loc[dorotheadb.target.isin(annTable.symbol)]
annTable = annTable.loc[annTable.symbol.isin(dorotheadb.target)]

# Prepare params

tfs = list(set(dorotheadb.loc[dorotheadb.tf.isin(tfs)].tf.tolist()))

print("List of TFs to use:", " ".join(tfs))

if len(tfs) == 0:
    print("Error")
    print("There are any TFs in the database")
    sys.exit()

universeGenes = list(set(annTable.symbol.tolist()))
targetGenes = list(set(dorotheadb.loc[dorotheadb.tf.isin(tfs)].target.tolist()))

biasFile = "data/"+database+"_"+annotation+".json"

if os.path.exists(biasFile):
    f = open(biasFile)
    bias = json.load(f)
    f.close()
else:
    bias = {gene:len(dorotheadb.loc[dorotheadb.target.isin([gene])].tf.tolist()) for gene in universeGenes}
    f = open(biasFile, "w")
    json.dump(bias, f)
    f.close()

de = {gene:1 if gene in targetGenes else 0 for gene in universeGenes}

dd = defaultdict(list)

for d in (bias, de):
    for key, value in d.items():
        dd[key].append(value)

bias = [dd[gene][0] for gene in universeGenes]
de = [dd[gene][1] for gene in universeGenes]

weigth = getGAModds(bias,de)
weigth = dict(zip(universeGenes, weigth))

annTableDict = annTable.groupby('annotation_id').apply(lambda x:x["symbol"].tolist()).to_dict()

terms = list(set(annTable.annotation_id))

N = len(list(set(annTable.symbol.tolist())))

print("Testing "+str(len(targetGenes))+" target genes in "+str(len(terms))+" terms.")
pvals = [calculate_pvalues(term,weigth,annTableDict,universeGenes,targetGenes,N) for term in terms]
pvals = pandas.concat(pvals)
pvals.to_csv(out,sep="\t",index=None)
