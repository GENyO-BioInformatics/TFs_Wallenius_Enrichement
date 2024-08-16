"
This script generates the common universe of Collectri with each annotation db.
This universe is understood as the set of TFs that at least each one targets one 
of the genes in the annotation db. See script figure_collectriAnnotationsCommonUniverse
and their resulting VennDiagrams in ./figures_data_n_scripts
"

library(glue)

sizes <- c(3, 10, 15, 20); size <- sizes[1]
set.seed(9)

outdir <- 'data/dbs_universes'
dir.create(outdir,showWarnings = F)

annotations <- c("KEGG","GO_BP","Reactome","WikiPathways"); annotation <- annotations[1]

organism <- 9606
collectri_raw_file <- glue("data/collectri_{organism}.tsv")
if (!file.exists(collectri_raw_file)){
  system("Rscript download_data.R")
}
collectriTFsGRN <- read.delim(collectri_raw_file,header = T,sep = '\t')

for (annotation in annotations){
  annotationFile <- glue("data/{annotation}.tsv")
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == organism,]
  
  collectriNannotDB <- glue("data/dbs_universes/TFs_universe_{annotation}_{organism}.txt")
  if (!file.exists(collectriNannotDB)){
    tfs <- c("IRF1","IRF2")
    write.table(tfs, file = "tfs.txt", row.names = F, quote = F, col.names = F)
    system(glue("python TFsEnrichment.py -ann {annotation} -tfs tfs.txt -m tfs"))
  }
  
  collectriNannotDB_TFsUniverse <- unique(read.delim(collectriNannotDB, header = F)[,1])
  collectriNannotDB_TFsUniverse <- collectriNannotDB_TFsUniverse[collectriNannotDB_TFsUniverse %in% collectriTFsGRN$tf & collectriNannotDB_TFsUniverse %in% annotationDF$symbol]
  collectriNannotDB <- glue("data/dbs_universes/TFs-Targets_universe_{annotation}_{organism}.txt")
  collectriNannotDB_TFsTargetsUniverse <- unique(read.delim(collectriNannotDB, header = F)[,1])
  
  for (size in sizes){
    randomOutdir <- paste0('random_lists_analysis/TFsLists/',annotation,'/size',size)
    dir.create(randomOutdir,recursive = T,showWarnings = F)
    for (n in 1:1000){
      TFsAnnoted <- sample(collectriNannotDB_TFsUniverse,size=size)
      write.table(TFsAnnoted, file = paste0(randomOutdir,'/',n,'_TFsAnnoted.txt'),
                  sep = '\t', row.names = F, col.names = F,quote = F)
      
      TFsTargetsAnnoted <- sample(collectriNannotDB_TFsTargetsUniverse,size=size)
      write.table(TFsTargetsAnnoted, file = paste0(randomOutdir,'/',n,'_TFsTargetsAnnoted.txt'),
                  sep = '\t', row.names = F, col.names = F,quote = F)
    }
  }
}
