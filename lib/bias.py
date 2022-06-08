import pandas, sys, os, urllib.request, random, collections, functools, operator
import json, getpass, glob, numpy, time, math, psutil, statistics
from scipy.stats import nchypergeom_wallenius
from multiprocessing import Pool
from pygam import LogisticGAM, s

def getPopAnnot(term,ann_table):
    return(len(list(set(ann_table.loc[ann_table.annotation_id.isin([term])].symbol.tolist()))))

def getAnnGenes(term,ann_table):
    return(list(set(ann_table[ann_table.annotation_id.isin([term])].symbol.tolist())))

def getInputAnnot(term,target_genes,annGenes):
    return(len([e for e in target_genes if e in annGenes]))

def by_database(database,ann_table):
    symbols = list(set(ann_table.symbol.tolist()))
    database_file = database + ".tsv"
    ### Read database with pandas
    database_dict = {}
    with open(os.path.join("processed_tables",database_file)) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i].rstrip().split("\t")
            if i != 0:
                tf = line[0]
                target = line[1]
                if target not in database_dict:
                    database_dict[target] = []
                database_dict[target].append(tf)
    database_bias = {target:len(database_dict[target]) for target in database_dict}
    for gene in symbols:
        if gene not in database_bias:
            database_bias[gene] = 0
    return database_bias

def getGAModds(lens,DEs):
    #lens,DEs = bias.bias,(bias.genes != '0').astype(int)
    x = numpy.array(lens)
    y = numpy.array(DEs)
    ww = numpy.argsort(x) # ww
    size = math.ceil(len(y) / 10)
    low = sum(y[ww][0:size])
    hi = sum(y[ww][(len(y) - size):len(y)])
    if hi <= low:
        reflectionFactor = 10^10
        x = reflectionFactor - x
        newX = x
    else:
        newX = x
        x = numpy.insert(x,0,0)
        y = numpy.insert(y,0,0)
    x = numpy.array([[x1] for x1 in x])
    y = numpy.array([[y1] for y1 in y])
    lams = numpy.exp(numpy.random.rand(100, 1))
    gam = LogisticGAM(s(0,n_splines=6,spline_order=3)).gridsearch(x,y,lam=lams,progress=False)
    probs = gam.predict_proba(newX)
    odds = probs / (1 - probs)
    return(odds)

def read_ann_file(annfile,tax_id):
    ann_table = {}
    with open(annfile) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i].rstrip().split("\t")
            if i != 0:
                annotation_id = line[1]
                symbol = line[0]
                organism = line[2]
                if organism == tax_id:
                    if annotation_id not in ann_table:
                        ann_table[annotation_id] = []
                    ann_table[annotation_id].append(symbol)
    return(ann_table)

from collections import ChainMap


def run_wallenius(ann,kegg_tax,bias,target_genes):
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


def by_term_wallenius(term,ann_table,target_genes,annotedbias,alpha,biasuniverse):
    ann_table_term = ann_table[term]
    term_genes = len(ann_table_term)
    genes_found = [gene for gene in ann_table_term if gene in target_genes]
    genes_found = len(genes_found)
    annotbias = annotedbias * (biasuniverse - genes_found) / (alpha - genes_found * annotedbias)
    annotbias = 1 if term_genes == genes_found else annotbias
    annotbias = {term:annotbias}
    return (annotbias)

def do_wallenius(inputAnnot,popAnnot,inputSize,popSize,odds):
    """
    popSize = all the target genes of a genome or system - M
    popAnnot = all the target genes of GO:XXXXXX - n
    inputSize = the input target genes - N
    inputAnnot = the input taget genes in GO:XXXXXX - x
    odds = the odds ratio of GO:XXXXXX - w 
    """
    pval = nchypergeom_wallenius(popSize,popAnnot,inputSize,odds).pmf(inputAnnot)
    return(pval)


annotation = ["KEGG","GO_BP","Reactome"]
# annotation = ["KEGG"]
database_used = ["dorothea_a_hsa","dorothea_b_hsa","dorothea_c_hsa","dorothea_d_hsa","dorothea_e_hsa"]
# database_used = ["dorothea_a_hsa"]

### Read TFs (tfs.txt)


database_to_use = database_used[0]

for database_to_use in database_used:
    print("Interactome "+database_to_use+" will be used")
    database_file = database_to_use + ".tsv"
    ### Read database with pandas
    database_pandas = pandas.read_csv(os.path.join("processed_tables",database_file),sep = "\t")
    organism = database_to_use.split("_")[-1]
    kegg_tax = tax_info[tax_info.kegg_id == organism]
    target_genes = list(set(database_pandas[database_pandas.tf.isin(tfs)].target.tolist()))
    ann = annotation[1]
    # bias = by_database(database_to_use)
    for ann in annotation:
        print("Doing "+ann)
        ann_table = pandas.read_csv("annotation/"+ann+".tsv",sep="\t")
        ann_table = ann_table[ann_table.organism == kegg_tax.taxonomy_id.tolist()[0]]
        print("Running empirical")
        # empirical_pvalue = run_empirical(ann,database_to_use)
        print("Running hypergeom")
        hypergeom_dataframe = run_hyper(ann_table,target_genes)
        print("Running wallenius")
        bias = by_database(database_to_use,ann_table)
        #wallenius_dataframe = run_wallenius(ann,kegg_tax,bias,target_genes)
        # empirical_pvalue.to_csv(database_to_use+"."+ann+".empirical_pvalue.csv",index=None)
        hypergeom_dataframe.to_csv(database_to_use+"."+ann+".hypergeom_pvalue.csv",index=None)
        wallenius_dataframe.to_csv(database_to_use+"."+ann+".wallenius_pvalue.csv",index=None)
        print("Done "+ann)
    print("Done Interactome "+database_to_use)
