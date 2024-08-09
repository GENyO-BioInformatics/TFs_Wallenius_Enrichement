try: 
  import pandas
  import os
  from collections import defaultdict
  import sys
  import argparse
  import numpy
  import pathlib
  from lib.calculate_stats import *
  from scipy.stats import nchypergeom_wallenius, hypergeom, rankdata
  import json
  from collections import Counter
  
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument("-tfs", "--tfs", type=str, help="list of tfs divided by commas (STAT1,STAT2,IRF1) or file to read in txt format",default = "STAT1,STAT2,IRF1")
  #parser.add_argument("-ann", "--annotation", type=str, help="select annotation to use", default="WikiPathways") # annotations = ["KEGG","GO_BP","Reactome","WikiPathways"]<
  #parser.add_argument("-out", "--output", type=str, help="output file", default="out.tsv")
  #parser.add_argument("-u", "--universe", type=str, help="txt file with the genes to define the background set")
  args = parser.parse_args()
  config = vars(args)
  
  tfsFile = config["tfs"]
  #annotation = config["annotation"]; #print("Perfoming enrichment analysis on "+annotation)
  #out = config["output"]
  #universe = config["universe"]
  
  # ##### FOR INTERACTIVE PROGRAMATIC EXAMPLE
  # tfsFile = 'STAT1,STAT2,IRF1'
  # tfsFile = 'data/TFsLists/WikiPathways/size3/169_TFsAnnoted.txt'
  organism = 9606
  annotation = tfsFile.split('/')[2]
  
  if tfsFile.endswith(".txt"):
    if os.path.exists(tfsFile):
        tfs = open(tfsFile).read().splitlines() 
    elif isinstance(tfsFile, str):
        tfs = tfsFile.split(",")
    elif ~isinstance(tfs, list):
      with open('errors_null_simulations.log','a') as loghandl:
        errorwhatever = "Incorrect format file. One column with TFs names is required"
        loghandl.write("ERROR {}:\n{}\n".format(tfsFile,errorwhatever))
        sys.exit()
  else:
    with open('errors_null_simulations.log','a') as loghandl:
      errorwhatever = "Incorrect input."
      loghandl.write("ERROR {}:\n{}\n".format(tfsFile,errorwhatever))
      sys.exit()
  
  print(tfsFile)
  #out = tfsFile.replace('TFsLists','EnrResults').replace('_TFsAnnoted.txt','.tsv')
  out = tfsFile.replace('TFsLists','EnrResults').replace('_TFsAnnoted.txt','_TFsTargetUniv.tsv')
  print(out)
  outdir = os.path.dirname(out)
  os.makedirs(outdir,exist_ok=True)
  
  # Read databases
  annFileTable = "data/"+annotation+".tsv"
  annTable = pandas.read_csv(annFileTable, sep="\t")
  annTable = annTable.loc[annTable.organism == organism]
  allTerms = list(set(annTable.annotation_id))
  
  collectridbfile = '/home/adriangarciamoreno/Desktop/TFs_Wallenius_Enrichement/data/collectri.tsv'
  collectridb = pandas.read_csv(collectridbfile,sep='\t')
  
  ## Stablish TFs UNIVERSES
  # PERFORM THE SEA ONLY USING THE TFs ASSOCIATED ANNOTATION UNIVERSE
  # this produced too low pvalues... 
  TFsUniverse = 'data/dbs_universes/{}-TFs_universe.txt'.format(annotation)
  TFsUniverse = open(TFsUniverse).read().splitlines() 
  # PERFORM THE SEA ONLY USING THE TFs-TARGETS ASSOCIATED ANNOTATION UNIVERSE
  TargetsTFsUniverse =  'data/dbs_universes/{}-TargetsTFs_universe.txt'.format(annotation)
  TargetsTFsUniverse = open(TargetsTFsUniverse).read().splitlines()
  # JUST IN CASE A TF is NOT INCLUDED IN TFs-TARGETS ASSOCIATED ANNOTATION UNIVERSE
  TFsUniverse = list(set(TFsUniverse + TargetsTFsUniverse))
  annTable_TFsUniverse =  annTable.loc[annTable.symbol.isin(TFsUniverse)]
  
  if len(set(TFsUniverse)) != len(set(annTable_TFsUniverse.symbol)):
    with open('errors_null_simulations.log','a') as loghandl:
      errorwhatever = "ERROR SETTING THE UNIVERSE TFs DIRECT"
      loghandl.write("ERROR {}:\n{}\n".format(tfsFile,errorwhatever))
      sys.exit()

  # TFsuniverseGenes == All the TFs from Collectri annotated in funtional DB
  TFsannTableDict = annTable_TFsUniverse.groupby('annotation_id').apply(lambda x:x["symbol"].tolist(), include_groups=False).to_dict()
  TFsN = len(list(set(annTable_TFsUniverse.symbol.tolist())))
  TFsterms = list(set(annTable_TFsUniverse.annotation_id))
  TFstotalTerms = len(TFsterms)
  TFsms = annTable_TFsUniverse.groupby('annotation_id').apply(lambda x:len(x["symbol"].tolist()), include_groups=False).to_dict()
  
  # N: universe || k: input size || xs: input genes in annotation || ms: all genes in annotation
  ### Analysis of TFs as genes || typeOfAnalysis == "TFs_Hypergeom"
  TFsk = len(tfs)
  if TFsk == 0:
      # res_tfs = dict(zip(terms,[[1,",".join(tfsTest)]]*totalTerms))
      pvals_TFs = [1]*totalTerms
  else:
      TFsxs = annTable_TFsUniverse.groupby('annotation_id').apply(lambda x: len(list(set(x["symbol"].tolist()) & set(tfs))), include_groups=False).to_dict() 
      if TFsms.keys() != TFsxs.keys():
        with open('errors_null_simulations.log','a') as loghandl:
          errorwhatever = 'UN-EXPECTED ERROR TFsms.keys() != TFsxs.keys()'
          loghandl.write("ERROR {}:\n{}\n".format(tfsFile,errorwhatever))
          sys.exit()

      # We substract 1 to x (set of input genes found in the annotation) 
      # because we are calculating the unliteral probabilty with 
      # the survival function that is defined as the probabily of finding at least > X-1 
      # genes in the annotation. Enrichment analyses done with the hpergeomentric 
      # distribution which is a discrete distribution. The alternative would be pmf(x) + sf(x). 
      # x=1; hypergeom.sf(x-1,1000,34,89) == hypergeom.pmf(x,1000,34,89) + hypergeom.sf(x,1000,34,89)
      # We directily use hypergeom.sf(x-1, ...) for code optimization reasons
      #arguments2test = list(zip(pandas.Series(TFsxs.values()) - 1, [TFsN]*TFstotalTerms, TFsms.values(), [TFsk]*TFstotalTerms))
      pvals_TFs = {key:hypergeom.sf(TFsxs[key] - 1,TFsN,TFsms[key],TFsk) for key,value in TFsxs.items()}
  
  #############################################################################################
  ### Analysis of TFs Target genes 
  ## GET TFs Targets UNIVERSE
  # tfsTargetsfile = tfsFile.replace('TFsAnnoted','TFsTargetsAnnoted')
  # if tfsTargetsfile.endswith(".txt"):
  #   if os.path.exists(tfsTargetsfile):
  #       tfsTargets = open(tfsTargetsfile).read().splitlines() 
  #   else:
  #       tfsTargets = tfsTargetsfile.split(",")
  #   if type(tfsTargets) != list:
  #     with open('errors_null_simulations.log','a') as loghandl:
  #       errorwhatever = "Incorrect format file. One column with TFs names is required"
  #       loghandl.write("ERROR {}:\n{}\n".format(tfsTargetsfile,errorwhatever))
  #       sys.exit()
  # else:
  #   with open('errors_null_simulations.log','a') as loghandl:
  #     errorwhatever = "Incorrect input."
  #     loghandl.write("ERROR {}:\n{}\n".format(tfsTargetsfile,errorwhatever))
  #     sys.exit()
  # 
  # ### Get the TFs Targets DB Universe
  # TargetsTFsUniverse =  'data/dbs_universes/{}-TargetsTFs_universe.txt'.format(annotation)
  # TargetsTFsUniverse = open(TargetsTFsUniverse).read().splitlines() 
  # annTable_TargetsTFsUniverse =  annTable.loc[annTable.symbol.isin(TargetsTFsUniverse)]
  # ### COLLECTRI MUST BE FILTERED FOR THE WALLENIOUS WEIGHT CALCULATION
  # subCollectri = collectridb.loc[collectridb.target.isin(TargetsTFsUniverse)]
  # 
  # if len(set(subCollectri.target)) != len(set(annTable_TargetsTFsUniverse.symbol)):
  #   with open('errors_null_simulations.log','a') as loghandl:
  #     errorwhatever = 'ERROR SETTING THE UNIVERSE TFs TARGETS'
  #     loghandl.write("ERROR {}:\n{}\n".format(tfsFile,errorwhatever))
  #     sys.exit()
  # 
  # ## Get TARGETS 
  # targetTest = list(set(subCollectri.loc[collectridb.tf.isin(tfsTargets)].target))
  # 
  # TargetsTFs_annTableDict = annTable_TargetsTFsUniverse.groupby('annotation_id').apply(lambda x:x["symbol"].tolist(), include_groups=False).to_dict()
  # TargetsTFs_N = len(list(set(annTable_TargetsTFsUniverse.symbol.tolist()))) # Universe
  # TargetsTFs_terms = list(set(annTable_TargetsTFsUniverse.annotation_id))
  # TargetsTFs_totalTerms = len(TargetsTFs_terms)
  # TargetsTFs_ms = annTable_TargetsTFsUniverse.groupby('annotation_id').apply(lambda x:len(x["symbol"].tolist()), include_groups=False).to_dict()
  # 
  # TargetsTFs_k = len(targetTest)
  # if TargetsTFs_k == 0:
  #     # NO TARGETS, NO HYPERGEOMETRIC NOR WALLENIUS CAN BE DONE
  #     pvals_TFsTargets = pvals_TFsWalle = [1]*totalTerms
  #     TargetsTFs_xs = dict(zip(annTable_TargetsTFsUniverse.groupby('annotation_id').groups.keys(), ['no_targets']*totalTerms)) 
  # else:
  #     TargetsTFs_xs = annTable_TargetsTFsUniverse.groupby('annotation_id').apply(lambda x:len(list(set(x["symbol"].tolist()) & set(targetTest))), include_groups=False).to_dict()
  #     if TargetsTFs_ms.keys() != TargetsTFs_xs.keys():
  #       with open('errors_null_simulations.log','a') as loghandl:
  #         errorwhatever = 'UN-EXPECTED ERROR TargetsTFs_ms.keys() != TargetsTFs_xs.keys()'
  #         loghandl.write("ERROR {}:\n{}\n".format(tfsFile,errorwhatever))
  #         sys.exit()
  #     
  #     ## typeOfAnalysis == "TFs_Hypergeom"    
  #     pvals_TFsTargets = {key:hypergeom.sf(TargetsTFs_xs[key] - 1,TargetsTFs_N,TargetsTFs_ms[key],TargetsTFs_k) 
  #                                          for key,value in TargetsTFs_xs.items()}
  # 
  #     ### "" || typeOfAnalysis == "Target_NonCentral"
  #     bias = Counter(subCollectri.target)
  #     universeGenesDF = pandas.DataFrame(index=list(set(annTable_TargetsTFsUniverse.symbol.tolist())))
  #     universeGenesDF['de'] = 0
  #     universeGenesDF.loc[targetTest,'de'] = 1
  #     universeGenesDF['bias'] = pandas.Series(bias)
  # 
  #     bias, de = universeGenesDF['bias'].to_list(), universeGenesDF['de'].to_list()
  #     universeGenesDF['weigth'] = getGAModds(bias, de)
  # 
  #     oddsRatios = {term:(numpy.mean(universeGenesDF.loc[TargetsTFs_annTableDict[term]].weigth) /
  #                   numpy.mean(universeGenesDF.loc[~universeGenesDF.index.isin(TargetsTFs_annTableDict[term])].weigth)) 
  #                   for term in TargetsTFs_xs.keys()}
  # 
  #     pvals_TFsWalle = {key:nchypergeom_wallenius.sf(TargetsTFs_xs[key] - 1,TargetsTFs_N,TargetsTFs_ms[key],TargetsTFs_k, oddsRatios[key])
  #                                                      for key,value in TargetsTFs_xs.items()}
  
  resultsDF = pandas.DataFrame({'terms':dict(zip(allTerms,allTerms)),
                                              'TFs_InputSize':TFsk, 
                                              'TFs_InputInAnnot':TFsxs, 
                                              'TFs_GenesAnnoted':TFsms, 
                                              'TFs_Universe':str(TFsN), 
                                              'TFs_Hypergeom':pvals_TFs},
                                              index=allTerms)
                                              # 'Targets_InputSize':TargetsTFs_k, 
                                              # 'Targets_InputInAnnot':TargetsTFs_xs, 
                                              # 'Targets_GenesAnnoted':TargetsTFs_ms, 
                                              # 'Targets_Universe':str(TargetsTFs_N), 
                                              # 'Target_Hypergeom':pvals_TFsTargets, 
                                              # 'oddsRatios':oddsRatios,
                                              # 'Target_NonCentral':pvals_TFsWalle,
                                              # 'Annotation':annotation},
  
  resultsDF.to_csv(out, sep="\t", index=None, header = True)
  
  with open('null_simulations.log','a') as loghandl:
      loghandl.write("DONE {}\n".format(tfsFile))

except Exception as errorwhatever:
    with open('errors_null_simulations.log','a') as loghandl:
      loghandl.write("ERROR {}:\n{}\n".format(tfsFile,errorwhatever))


# .sf(TFs_InputInAnnot - 1, TFs_Universe, TFs_InputSize, TFs_GenesAnnoted)
# hypergeom.sf(2-1,680,3,24)
# hypergeom.sf(12-1,4302,132,175)
# nchypergeom_wallenius.sf(12-1,4302,132,175,1.71455454373074)



