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

# Establish the size length of TFs
tfSize = [2,3,4,5,10,20,30]

# Using 1000 repetitions take a long time
repeats = 1000

# Select the annotation to use
annotations = ["KEGG","GO_BP","Reactome","WikiPathways"]
dorotheadb = getDorotheaDB(confidence="C", org=9606)

def multiprocessingRepeats(params):
    repeat, randomTFs, annotation, size = params
    tfsUse = ",".join(randomTFs)
    outFile = "results/"+annotation+"_"+str(size)+"_"+str(repeat + 1)+".tsv"
    os.system("python3 TFsWalleniusEnrichement.py --tfs "+tfsUse+" --output "+outFile+" --annotation "+annotation)

annotation = "KEGG"
t1_start = process_time()
for annotation in annotations:
    print(annotation)
    annFileTable = "data/"+annotation+".tsv"
    annTable = pandas.read_csv(annFileTable, sep="\t")
    annTable = annTable.loc[annTable.organism == 9606]

    # Filter databases
    dorotheadbAnn = dorotheadb.loc[dorotheadb.target.isin(annTable.symbol)]
    annTable = annTable.loc[annTable.symbol.isin(dorotheadbAnn.target)]
    tfs = list(set(dorotheadbAnn.tf.tolist()))

    universeGenes = list(set(annTable.symbol.tolist()))
    biasFile = "data/C_"+annotation+".json"

    if os.path.exists(biasFile):
        f = open(biasFile)
        bias = json.load(f)
        f.close()
    else:
        bias = {gene:len(dorotheadbAnn.loc[dorotheadbAnn.target.isin([gene])].tf.tolist()) for gene in universeGenes}
        f = open(biasFile, "w")
        json.dump(bias, f)
        f.close()

    for size in tfSize:
        random.seed(size)
        randomTFs = [random.sample(tfs,size) for i in range(repeats)]
        input_list = [(repeat, randomTFs[repeat], annotation, size) for repeat in range(repeats)]
        with multiprocessing.Pool() as pool:
            catch_res = pool.map(multiprocessingRepeats,input_list)

        resultFiles = glob.glob("results/"+annotation+"_"+str(size)+"_*")
        content = pandas.concat([pandas.read_csv(fileName, sep = "\t") for fileName in resultFiles])
        catch_res = [os.remove(fileName) for fileName in resultFiles]
        terms = list(set(content.term.tolist()))
        check_hyper = {term:numpy.count_nonzero(content[content.term == term].pval_hyper.to_numpy() < 0.05) for term in terms}
        check_wall = {term:numpy.count_nonzero(content[content.term == term].pval_wall.to_numpy() < 0.05) for term in terms}

        dd = defaultdict(list)
        for d in (check_hyper, check_wall):
            for key, value in d.items():
                dd[key].append(value / repeats)

        n_dict = [pandas.DataFrame({"term":[term],"hypergeom":[dd[term][0]],"wallenius":[dd[term][1]]}) for term in terms]
        results = pandas.concat(n_dict)

        outFileName = "results/"+annotation+"_"+str(size)+".tsv"
        results.to_csv(outFileName, sep = "\t", index = None)

t1_stop = process_time()
print("Elapsed time during the whole program in seconds:",t1_stop-t1_start)
