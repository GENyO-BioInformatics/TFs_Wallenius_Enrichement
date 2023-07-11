from lib.dorotheadb import *
import random
import os
import glob
import numpy
import json
import multiprocessing
from scipy.stats import combine_pvalues, rankdata
import pandas

# Establish the size length of TFs
tfSize = [2,5,10,20,30]

# Using 1000 repetitions take a long time
repeats = 1000

# Select the annotation to use
annotations = ["KEGG","GO_BP","Reactome","WikiPathways"]
dorotheadb = getDorotheaDB(confidence="C", org=9606)
tfs = list(set(dorotheadb.tf.tolist()))

def multiprocessingRepeats(params):
    repeat, randomTFs, annotation, size = params
    tfsUse = ",".join(randomTFs)
    outFile = "null_simulations/"+annotation+"_"+str(size)+"_"+str(repeat + 1)+".json"
    os.system("python3 TFsEnrichment.py --tfs "+tfsUse+" --output "+outFile+" --annotation "+annotation)

for annotation in annotations:
    print(annotation)

    for size in tfSize:
        random.seed(size)
        randomTFs = [random.sample(tfs, size) for i in range(repeats)]
        input_list = [(repeat, randomTFs[repeat], annotation, size) for repeat in range(repeats)]
        # with multiprocessing.Pool() as pool:
        #     catch_res = pool.map(multiprocessingRepeats,input_list)

        resultFiles = glob.glob("null_simulations/"+annotation+"_"+str(size)+"_*.json")
        fileName = resultFiles[0]
        listOfResults = {}
        for fileName in resultFiles:
            f = open(fileName)
            results = json.load(f)
            f.close()

            typesOfAnalysis = results.keys()
            if len(listOfResults.keys()) == 0:
                listOfResults = {typeOfAnalysis:{"pval":{}, "w":{}} for typeOfAnalysis in typesOfAnalysis}
            for typeOfAnalysis in typesOfAnalysis:
                if len(listOfResults[typeOfAnalysis]["pval"]) == 0:
                    listOfResults[typeOfAnalysis]["pval"] = {term:[results[typeOfAnalysis]["pval"][term]] for term in results[typeOfAnalysis]["pval"].keys()}
                    listOfResults[typeOfAnalysis]["ranking"] = {term:[results[typeOfAnalysis]["ranking"][term]] for term in results[typeOfAnalysis]["pval"].keys()}

                else:
                    for term in results[typeOfAnalysis]["pval"].keys():
                        listOfResults[typeOfAnalysis]["pval"][term].append(results[typeOfAnalysis]["pval"][term])
                        listOfResults[typeOfAnalysis]["ranking"][term].append(results[typeOfAnalysis]["ranking"][term])

        # catch_res = [os.remove(fileName) for fileName in resultFiles]
        stats = pandas.DataFrame({"annotation_id":[],"typeOfAnalysis":[],"p_times":[],"ranking":[]})
        for typeOfAnalysis in typesOfAnalysis:
            typeOfAnalysisRes = listOfResults[typeOfAnalysis]
            terms = list(typeOfAnalysisRes["pval"].keys())
            resComb = numpy.array([sum(numpy.array(typeOfAnalysisRes["pval"][term]) < 0.05) / repeats for term in terms])
            resRank = numpy.array([numpy.median(numpy.array(typeOfAnalysisRes["ranking"][term])) for term in terms])
            statsDF = pandas.DataFrame({"annotation_id":terms, "typeOfAnalysis":[typeOfAnalysis]*len(terms), "p_times":resComb, "ranking":resRank}).sort_values(by = "ranking")
            stats = pandas.concat([stats,statsDF])

        concatContent = stats.sort_values(by="ranking")

        outFileName = "random_list_analysis/"+annotation+"_"+str(size)+".tsv"
        concatContent.to_csv(outFileName,sep="\t",index=None)
