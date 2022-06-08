import pandas, os, subprocess, pyreadr
from collections import Counter
from lib.dorotheadb import *

dorotheadb = getDorotheaDB()

dorotheadb

dorotheadb

#annotation = ["KEGG","GO_BP","Reactome"]
annotation = "KEGG"
tfs = ['STAT1','STAT2','IRF1','IRF3']

tfsPERtarget = Counter(dorotheadb.target)

tfsPERtarget

    results = run_wallenius(ann,kegg_tax,bias,target_genes)




tax_id = str(kegg_tax.taxonomy_id.tolist()[0])
annfile = "annotation/"+ann+".tsv"
ann_table = read_ann_file(annfile,tax_id)

engene = {"id":[],"bias":[],"genes":[]}
for gene in bias:
    engene["id"].append(gene)
    engene["bias"].append(bias[gene])
    if gene in target_genes:
        engene["genes"].append(gene)
    else:
        engene["genes"].append("0")
engene = pandas.DataFrame(engene)
odds = getGAModds(engene.bias,(engene.genes != '0').astype(int))
engene['odds'] = odds
alpha = sum(engene.odds)
biasuniverse = len(engene.index)
annotedbias = numpy.mean(engene.loc[engene.genes != "0"].odds)

results = [by_term_wallenius(term,ann_table,target_genes,annotedbias,alpha,biasuniverse) for term in list(ann_table)]
results = dict(ChainMap(*results))

terms = list(ann_table)
popSize = len(engene.index)
inputSize = len(target_genes)

ann_table = pandas.read_csv("annotation/"+ann+".tsv",sep="\t")
ann_table = ann_table[ann_table.organism == kegg_tax.taxonomy_id.tolist()[0]]

with Pool() as pool:
    popAnnot = pool.starmap(getPopAnnot,zip(terms,[ann_table]*len(terms)))
with Pool() as pool:
    annGenes = pool.starmap(getAnnGenes,zip(terms,[ann_table]*len(terms)))
with Pool() as pool:
    inputAnnot = pool.starmap(getInputAnnot,zip(terms,[target_genes]*len(terms),annGenes))
i = list(range(len(terms)))
annotbias = [results[term] for term in terms]
input_list = [(inputAnnot[e],popAnnot[e],inputSize,popSize,annotbias[e]) for e in i]
with Pool() as pool:
    pvalues = pool.starmap(do_wallenius,input_list)

wallenius_dataframe = pandas.DataFrame({"term":terms,"pvalue":pvalues})
wallenius_dataframe = wallenius_dataframe.sort_values(by=["pvalue"])
return(wallenius_dataframe)
