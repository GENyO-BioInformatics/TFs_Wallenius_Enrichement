import numpy
import math
from pygam import LogisticGAM, s
from scipy.stats import nchypergeom_wallenius, hypergeom
import pandas

def getGAModds(lens,DEs):
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
    gam = LogisticGAM(s(0,n_splines=6,spline_order=3,constraints='monotonic_inc')).fit(x,y)
    probs = gam.predict_proba(newX)
    return(probs)

def calculate_pvalues(term,weigth,annTableDict,universeGenes,targetGenes,N):
    insideGenes = annTableDict[term]
    insideWeigthSum = numpy.mean([weigth[gene] for gene in insideGenes])
    outsideGenes = [gene for gene in universeGenes if gene not in insideGenes]
    outsideWeigthSum = numpy.mean([weigth[gene] for gene in outsideGenes])
    odd_ratio = insideWeigthSum / outsideWeigthSum

    k = len(targetGenes)
    m = len(insideGenes)
    x = len(list(set(insideGenes) & set(targetGenes)))

    pval_hyper = hypergeom.sf(x-1,N,m,k)
    if odd_ratio == 1:
        pval_wall = pval_hyper
    else:
        pval_wall = nchypergeom_wallenius.sf(x-1, N, m, k, odd_ratio)

    if pval_wall == 0:
        pval_wall = pval_hyper
    res = pandas.DataFrame({"term":[term],"pval_hyper":[pval_hyper],"pval_wall":[pval_wall]})
    return(res)
