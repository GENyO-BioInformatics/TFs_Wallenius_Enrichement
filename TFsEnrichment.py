import time
import pandas as pd
import numpy as np
import argparse
import os
import sys
from scipy.stats import nchypergeom_wallenius, hypergeom, rankdata
from collections import Counter
from lib.calculate_stats import getGAModds
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# python TFsEnrichment.py -ann KEGG -tfs tfs.txt -out out.txt -m targetsW
parser.add_argument("-org", "--organism", type=int, help="select tax id: 9606 for human and 10090 for mouse", default=9606)
parser.add_argument("-ann", "--annotation", type=str, help="select annotation to use", default="GO_BP")
parser.add_argument("-out", "--output", type=str, help="output file", default="out.txt")
parser.add_argument("-tfs", "--tfs", type=str, help="text file with TFs", required=True)
parser.add_argument("-m", "--method", type=str, help="select method to use: tfs, targetsF, targetsW", default="targetsW")
args = parser.parse_args()
config = vars(args)

tfsFile = config["tfs"]
out = config["output"]
annotation = config["annotation"]
organism = config["organism"]
method = config["method"]

# Check if tfsFile exists and has the right format
try:
    with open(tfsFile, 'r') as file:
        tfs = file.read().splitlines()
        if len(tfs) == 0:
            raise ValueError(f"The file {os.path.basename(tfsFile)} is empty. It must have the following structure:\n\nSTAT1\nSTAT2\nIRF1\n")
except FileNotFoundError:
    print(f"\nERROR: The file '{os.path.basename(tfsFile)}' does not exist. Please ensure that the file is located at '{tfsFile}' and contains three lines with the following structure:\n\nSTAT1\nSTAT2\nIRF1\n")
except ValueError as e:
    print(f"\nERROR: {e}")

# Set unique tfs (in case there are duplicated entries)
tfs = list(set(tfs))

# Check if exists annotation file and the organism has annotations
annFileTable = f"data/{annotation}.tsv"
try:
    annTable = pd.read_csv(annFileTable, sep="\t")
    if organism not in annTable["organism"].values:
        raise ValueError(f"Organism {organism} does not have annotations in the file {os.path.basename(annFileTable)}. We only have annotations for 9606 (human) and 10090 (mouse)")
    annTable = annTable.loc[annTable.organism == organism]
except FileNotFoundError:
    print(f"\nERROR: The file '{os.path.basename(annFileTable)}' does not exist. Please ensure that the file is located at '{annFileTable}'. It must be in data folder (data/KEGG.tsv for example)\n")
except ValueError as e:
    print(f"\nERROR: {e}\n")

# check if exists collectri file
collectriDBFile = f"data/collectri_{organism}.tsv"

if not os.path.exists(collectriDBFile):
    print(f"Preparing database. It will take a while but it will only be done once. You can check the file at '{collectriDBFile}'")
    os.system(f"Rscript download_data.R -o {organism}")

collectriDB = pd.read_csv(collectriDBFile,sep='\t')

# check if universe-db files exists
tfsUniverse = f"data/dbs_universes/TFs_universe_{annotation}_{organism}.txt" # TF as symbol in annotation
targetsTFsUniverse = f"data/dbs_universes/Targets-TFs_universe_{annotation}_{organism}.txt" # target as symbol in annotation
tfsTargetsUniverse = f"data/dbs_universes/TFs-Targets_universe_{annotation}_{organism}.txt" # TFs of targets as symbol in annotation

def prepare_universe(organism, annTable, tfsUniverse, tfsTargetsUniverse, targetsTFsUniverse):
    os.makedirs("data/dbs_universes",exist_ok=True)
    collectri_raw = pd.read_csv(f"data/collectri_raw_{organism}.tsv",sep='\t')
    # get subset if collectri TFs annotated in annTable
    subcollectri = collectri_raw[collectri_raw.source_genesymbol.isin(annTable.symbol)]
    subcollectri_tfs = list(set(collectri_raw.source_genesymbol.to_list()))
    subcollectri_targets = list(set(collectri_raw.target_genesymbol.to_list()))
    tfsUniverseList = list(set(subcollectri_tfs + subcollectri_targets)) # universe for TFs are all TFs and targets
    # write to file
    with open(tfsUniverse,'w') as f:
        f.write('\n'.join(tfsUniverseList))
    # get subset if collectri targets annotated in annTable
    subcollectri = collectri_raw[collectri_raw.target_genesymbol.isin(annTable.symbol)]
    subcollectri_tfs = list(set(subcollectri.source_genesymbol.to_list())) # list of TFs that have a target
    # write to file
    with open(tfsTargetsUniverse,'w') as f:
        f.write('\n'.join(subcollectri_tfs))
    # get the set of targets
    subcollectri_targets = list(set(subcollectri.target_genesymbol.to_list())) # list of targets
    # write to file
    with open(targetsTFsUniverse,'w') as f:
        f.write('\n'.join(subcollectri_targets))

    

if not os.path.exists(tfsUniverse) or not os.path.exists(targetsTFsUniverse) or not os.path.exists(tfsTargetsUniverse):
    print(f"Preparing universe database. It will take a while but it will only be done once.")
    prepare_universe(organism, annTable, tfsUniverse, tfsTargetsUniverse, targetsTFsUniverse)

tfsUniverse = open(tfsUniverse).read().splitlines()
tfsTargetsUniverse = open(tfsTargetsUniverse).read().splitlines()
targetsTFsUniverse = open(targetsTFsUniverse).read().splitlines()

def tf_fisher(annTable, tfs, tfsUniverse):
    annTable_TFsUniverse =  annTable.loc[annTable.symbol.isin(tfsUniverse)]
    # N: universe || k: input size || xs: input genes in annotation || ms: all genes in annotation
    k = len(tfs)
    N = len(list(set(annTable_TFsUniverse.symbol.tolist())))
    xs = annTable_TFsUniverse.groupby('annotation_id').apply(lambda x: len(list(set(x["symbol"].tolist()) & set(tfs)))).to_dict()
    ms = annTable_TFsUniverse.groupby('annotation_id').apply(lambda x:len(x["symbol"].tolist())).to_dict()
    # We substract 1 to x (set of input genes found in the annotation) 
    # because we are calculating the unliteral probability with 
    # the survival function that is defined as the probability of finding at least > X-1 
    # genes in the annotation. Enrichment analyses done with the hpergeomentric 
    # distribution which is a discrete distribution. The alternative would be pmf(x) + sf(x). 
    # x=1; hypergeom.sf(x-1,1000,34,89) == hypergeom.pmf(x,1000,34,89) + hypergeom.sf(x,1000,34,89)
    # We directily use hypergeom.sf(x-1, ...) for code optimization reasons
    pvals_TFs = {key:hypergeom.sf(xs[key] - 1,N,ms[key],k) for key,value in xs.items()}
    return pvals_TFs

def target_fisher(annTable, tfs, targetsTFsUniverse, collectriDB):
    annTable_TargetsTFsUniverse =  annTable.loc[annTable.symbol.isin(targetsTFsUniverse)]
    # filter collectri
    subCollectri = collectriDB.loc[collectriDB.target.isin(targetsTFsUniverse)]
    ## Get TARGETS 
    targetTest = list(set(subCollectri.loc[collectriDB.tf.isin(tfs)].target))
    # N: universe || k: input size || xs: input genes in annotation || ms: all genes in annotation
    N = len(list(set(annTable_TargetsTFsUniverse.symbol.tolist())))
    k = len(targetTest)
    xs = annTable_TargetsTFsUniverse.groupby('annotation_id').apply(lambda x:len(list(set(x["symbol"].tolist()) & set(targetTest)))).to_dict()
    ms = annTable_TargetsTFsUniverse.groupby('annotation_id').apply(lambda x:len(x["symbol"].tolist())).to_dict()
    pvals_targets = {key:hypergeom.sf(xs[key] - 1,N,ms[key],k) 
                                           for key,value in xs.items()}
    return pvals_targets

def target_wallenius(annTable, tfs, targetsTFsUniverse, collectriDB):
    annTable_TargetsTFsUniverse =  annTable.loc[annTable.symbol.isin(targetsTFsUniverse)]
    targetsTFs_annTableDict = annTable_TargetsTFsUniverse.groupby('annotation_id').apply(lambda x:x["symbol"].tolist()).to_dict()
    # filter collectri
    subCollectri = collectriDB.loc[collectriDB.target.isin(targetsTFsUniverse)]
    targetTest = list(set(subCollectri.loc[collectriDB.tf.isin(tfs)].target))
    bias = Counter(subCollectri.target)
    universeGenesDF = pd.DataFrame(index=list(set(annTable_TargetsTFsUniverse.symbol.tolist())))
    universeGenesDF['de'] = 0
    universeGenesDF.loc[targetTest,'de'] = 1
    universeGenesDF['bias'] = pd.Series(bias)
    bias, de = universeGenesDF['bias'].to_list(), universeGenesDF['de'].to_list()
    universeGenesDF['weigth'] = getGAModds(bias, de)
    N = len(list(set(annTable_TargetsTFsUniverse.symbol.tolist())))
    k = len(targetTest)
    xs = annTable_TargetsTFsUniverse.groupby('annotation_id').apply(lambda x:len(list(set(x["symbol"].tolist()) & set(targetTest)))).to_dict()
    ms = annTable_TargetsTFsUniverse.groupby('annotation_id').apply(lambda x:len(x["symbol"].tolist())).to_dict()
    oddsRatios = {term:(np.mean(universeGenesDF.loc[targetsTFs_annTableDict[term]].weigth) /
                    np.mean(universeGenesDF.loc[~universeGenesDF.index.isin(targetsTFs_annTableDict[term])].weigth)) 
                    for term in xs.keys()}
    pvals_TFsWalle = {key:nchypergeom_wallenius.sf(xs[key] - 1,N,ms[key],k, oddsRatios[key])
                                                       for key,value in xs.items()}
    return pvals_TFsWalle

terms = list(set(annTable.annotation_id))

if method == "tfs":
    pvals_TFsFisher = tf_fisher(annTable, tfs, tfsUniverse)
    pvals_TFsFisher = [pvals_TFsFisher[key] if key in pvals_TFsFisher else 1 for key in terms]
    results = pd.DataFrame({'term':terms,'pvals_TFsFisher':pvals_TFsFisher, "annotation": annotation, "k": len(tfs)})
    results.to_csv(out, sep = "\t", index = None)

elif method == "targetsF":
    pvals_targetsFisher = target_fisher(annTable, tfs, targetsTFsUniverse, collectriDB)
    pvals_targetsFisher = [pvals_targetsFisher[key] if key in pvals_targetsFisher else 1 for key in terms]
    results = pd.DataFrame({'term':terms,'pvals_targetsFisher':pvals_targetsFisher, "annotation": annotation, "k": len(tfs)})
    results.to_csv(out, sep = "\t", index = None)

elif method == "targetsW":
  pvals_targetsWallenius = target_wallenius(annTable, tfs, targetsTFsUniverse, collectriDB)
  pvals_targetsWallenius = [pvals_targetsWallenius[key] if key in pvals_targetsWallenius else 1 for key in terms]
  results = pd.DataFrame({'term':terms,'pvals_targetsWallenius':pvals_targetsWallenius, "annotation": annotation, "k": len(tfs)})
  results.to_csv(out, sep = "\t", index = None)