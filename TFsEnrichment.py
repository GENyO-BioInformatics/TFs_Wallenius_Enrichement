import pandas
import os
from collections import defaultdict
import sys
import argparse
import numpy
import pathlib
from lib.dorotheadb import *
from lib.calculate_stats import *
from scipy.stats import nchypergeom_wallenius, hypergeom, rankdata
import json

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-db", "--database", help="select level of confidence for dorothea", default="C")
parser.add_argument("-org", "--organism", type=int, help="select tax id: 9606 for human and 10090 for mouse", default=9606)
parser.add_argument("-ann", "--annotation", type=str, help="select annotation to use", default="WikiPathways")
parser.add_argument("-out", "--output", type=str, help="output file", default="out.json")
parser.add_argument("-tfs", "--tfs", type=str, help="list of tfs divided by commas (STAT1,STAT2,IRF1) or file to read in txt format",default = "STAT1,STAT2,IRF1")
args = parser.parse_args()
config = vars(args)

tfs = config["tfs"]
out = config["output"]
file_extension = pathlib.Path(out).suffix
annotation = config["annotation"]
database = config["database"]
organism = config["organism"]

# print("Perfoming enrichment analysis on "+annotation)

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

## Change the next line if you want to perform less types of enrichment analysis
typesOfAnalysis = ["TFs_Hypergeom", "Target_Hypergeom", "Target_NonCentral"]
universeGenes = list(set(annTable.symbol.tolist()))
annTableDict = annTable.groupby('annotation_id').apply(lambda x:x["symbol"].tolist()).to_dict()
terms = list(set(annTable.annotation_id))
N = len(list(set(annTable.symbol.tolist())))



if file_extension == ".json":
    fullRes = {}
else:
    fullRes = []

for typeOfAnalysis in typesOfAnalysis:
    # print("Performing "+typeOfAnalysis+" enrichment analysis")
    if typeOfAnalysis == "TFs_Hypergeom":
        dorotheadbTest = dorotheadb.loc[dorotheadb.tf.isin(annTable.symbol)]
        tfsTest = list(set(dorotheadbTest.loc[dorotheadbTest.tf.isin(tfs)].tf.tolist()))
        res = {}
        if len(tfsTest) == 0:
            for term in terms:
                res[term] = [1,",".join(tfsTest)]
        else:
            for term in terms:
                insideGenes = annTableDict[term]
                k = len(tfsTest)
                m = len(insideGenes)
                x = len(list(set(insideGenes) & set(tfsTest)))
                pval = hypergeom.sf(x-1,N,m,k)
                res[term] = [pval,",".join(tfsTest)]
    elif typeOfAnalysis == "Target_Hypergeom":
        dorotheadbTest = dorotheadb.loc[dorotheadb.target.isin(annTable.symbol)]
        targetTest = list(set(dorotheadbTest.loc[dorotheadbTest.tf.isin(tfs)].target.tolist()))
        res = {}
        if len(targetTest) == 0:
            for term in terms:
                res[term] = [1,",".join(targetTest)]
        else:
            for term in terms:
                insideGenes = annTableDict[term]
                k = len(targetTest)
                m = len(insideGenes)
                x = len(list(set(insideGenes) & set(targetTest)))
                pval = hypergeom.sf(x-1,N,m,k)
                res[term] = [pval,",".join(targetTest)]
    else:
        dorotheadbTest = dorotheadb.loc[dorotheadb.target.isin(annTable.symbol)]
        targetTest = list(set(dorotheadbTest.loc[dorotheadbTest.tf.isin(tfs)].target.tolist()))
        res = {}
        if len(targetTest) == 0:
            for term in terms:
                res[term] = [1,",".join(targetTest)]
        else:
            biasFile = "data/"+database+"_"+annotation+".json"
            if os.path.exists(biasFile):
                f = open(biasFile)
                bias = json.load(f)
                f.close()
            else:
                bias = {gene:len(dorotheadbTest.loc[dorotheadbTest.target.isin([gene])].tf.tolist()) for gene in universeGenes}
                f = open(biasFile, "w")
                json.dump(bias, f)
                f.close()
            de = {gene:1 if gene in targetTest else 0 for gene in universeGenes}
            dd = defaultdict(list)
            for d in (bias, de):
                for key, value in d.items():
                    dd[key].append(value)

            bias = [dd[gene][0] for gene in universeGenes]
            de = [dd[gene][1] for gene in universeGenes]
            weigth = getGAModds(bias, de)
            weigth = dict(zip(universeGenes, weigth))

            for term in terms:
                insideGenes = annTableDict[term]
                insideWeigthSum = numpy.mean([weigth[gene] for gene in insideGenes])
                outsideGenes = [gene for gene in universeGenes if gene not in insideGenes]
                outsideWeigthSum = numpy.mean([weigth[gene] for gene in outsideGenes])
                oddRatio = insideWeigthSum / outsideWeigthSum
                k = len(targetTest)
                m = len(insideGenes)
                x = len(list(set(insideGenes) & set(targetTest)))
                if oddRatio == 1:
                    pval = hypergeom.sf(x-1,N,m,k)
                else:
                    pval = nchypergeom_wallenius.sf(x-1, N, m, k, oddRatio)
                    if pval == 0:
                        pval = hypergeom.sf(x-1,N,m,k)
                res[term] = [pval,",".join(targetTest)]
    pvals = [res[term][0] for term in terms]
    input = [res[term][1] for term in terms]
    ranking = rankdata(pvals)
    ranking = {terms[i]: ranking[i] for i in range(len(terms))}
    pvals = {terms[i]: pvals[i] for i in range(len(terms))}
    input = {terms[i]: input[i] for i in range(len(terms))}
    if file_extension == ".json":
        fullRes[typeOfAnalysis] = {"pval":pvals, "ranking":ranking, "input": input}
    else:
        fullRes.append({"annotation_id":terms, "typeOfAnalysis":[typeOfAnalysis]*len(terms), "pval":list(pvals.values()), "ranking":list(ranking.values()), "input":list(input.values())})

if file_extension == ".json":
    f = open(out, "w")
    json.dump(fullRes, f)
    f.close()
else:
    concatContent = pandas.DataFrame({"annotation_id":[],"typeOfAnalysis":[],"pval":[],"ranking":[],"input":[]})
    for i in range(len(fullRes)):
        content = pandas.DataFrame.from_dict(fullRes[i])
        concatContent = pandas.concat([concatContent,content])
    concatContent = concatContent.sort_values(by="ranking")
    if (file_extension == ".txt") | (file_extension == ".tsv"):
        concatContent.to_csv(out,sep="\t",index=None)
    else:
        concatContent.to_csv(out,sep=",",index=None)
