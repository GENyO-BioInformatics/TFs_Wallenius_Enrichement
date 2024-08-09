"
This script generates the common universe of Collectri with each annotation db.
This universe is understood as the set of TFs that at least each one targets one 
of the genes in the annotation db. See script figure_collectriAnnotationsCommonUniverse
and their resulting VennDiagrams in ./figures_data_n_scripts
"
sizes <- c(3, 10, 15, 20); size <- sizes[1]
set.seed(9)

outdir <- 'data/dbs_universes'
dir.create(outdir,showWarnings = F)

collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
statisticsDF <- c("BaseUniverse","Collectri",length(unique(collectriTFsGRN$source)),
                  length(unique(collectriTFsGRN$target)),nrow(collectriTFsGRN))
annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv')
annotationFile <- annotationFiles[1]

for (annotationFile in annotationFiles){
  annotation <- gsub('.tsv','',basename(annotationFile))
  print(annotationFile)
  
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == 9606,]
  
  ### GET SUBSET OF COLLECTRI TFs THAT ARE ANNOTATED IN DB,
  ### THE SET OF TFs FROM HERE MUST BE USED TO LIMIT THE ANNOTATION DBs
  ### WHEN TESTING DIRECTLY THE TFs ANNOTATIONS
  subCollectri <- collectriTFsGRN[collectriTFsGRN$source_genesymbol %in% annotationDF$symbol,] 
  statisticsDF <- rbind(statisticsDF,cbind("TFsUniverse",annotation,length(unique(subCollectri$source_genesymbol)),
                                           length(unique(subCollectri$target_genesymbol)),nrow(subCollectri)))
  collectriNannotDB_TFsUniverse <- unique(subCollectri$source_genesymbol)
  
  TFsOutfile <- gsub('.tsv','-TFs_universe.txt',gsub('data','data/dbs_universes',annotationFile))
  write(collectriNannotDB_TFsUniverse,TFsOutfile,sep = '\n')
  
  
  ### GET SUBSET OF COLLECTRI WHOSE TARGETS THAT ARE IN THE ANNOTATION DB
  subCollectri <- collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% annotationDF$symbol,] 
  length(unique(subCollectri$target_genesymbol))
  statisticsDF <- rbind(statisticsDF,cbind("TargetsUniverse",annotation,length(unique(subCollectri$source_genesymbol)),
                                           length(unique(subCollectri$target_genesymbol)),nrow(subCollectri)))
  
  ### THE SET OF TFs FROM HERE MUST BE USED FOR THE RANDOM LISTS OF TFs TO
  ### PERFORM THE TEST BASED ON TARGETS
  collectriNannotDB_TFsTargetsUniverse <- unique(subCollectri$source_genesymbol)
  annotation <- gsub('.tsv','',basename(annotationFile))
  TFsTargetsUniverseOutfile <- gsub('.tsv','-TFsTargets_universe.txt',gsub('data','data/dbs_universes',annotationFile))
  write(collectriNannotDB_TFsTargetsUniverse,TFsTargetsUniverseOutfile,sep = '\n')

  ### THE SET OF COLLECTRI TARGETS FROM HERE MUST BE USED TO LIMIT THE ANNOTATION DBs
  ### AND THEN THE COLLECTRI UNIVERSE WHEN TESTING VIA TARGET GENES
  collectriNannotDB_TargetsTFsUniverse <- unique(subCollectri$target_genesymbol)
  annotation <- gsub('.tsv','',basename(annotationFile))
  TargetsTFsUniverseOutfile <- gsub('.tsv','-TargetsTFs_universe.txt',gsub('data','data/dbs_universes',annotationFile))
  write(collectriNannotDB_TargetsTFsUniverse,TargetsTFsUniverseOutfile,sep = '\n')

  for (size in sizes){
    randomOutdir <- paste0('data/TFsLists/',annotation,'/size',size)
    dir.create(randomOutdir,recursive = T,showWarnings = F)
    for (n in 1:1000){
      TFsAnnoted <- sample(collectriNannotDB_TFsUniverse,size=size,replace = T)
      write.table(TFsAnnoted, file = paste0(randomOutdir,'/',n,'_TFsAnnoted.txt'),
                  sep = '\t', row.names = F, col.names = F,quote = F)

      TFsTargetsAnnoted <- sample(collectriNannotDB_TFsTargetsUniverse,size=size,replace = T)
      write.table(TFsTargetsAnnoted, file = paste0(randomOutdir,'/',n,'_TFsTargetsAnnoted.txt'),
                  sep = '\t', row.names = F, col.names = F,quote = F)
    }
  }
}

statisticsDF <- as.data.frame(statisticsDF)
colnames(statisticsDF) <- c('UniverseLevel','Annnotation','TFs','Targets','Interactions')
write.table(statisticsDF,file.path(outdir,"statistics.tsv"),sep = '\t',row.names = F,col.names = T, quote = F)


### CURRENT ISSUES
## 1
# In Collectri some TFs are noted as Complexes
# This are not annotated as such in functional databases
# Since random TFs are picked from collectri some TFs will not be analysed
# in the first approach

