
TFregulons <- read.delim("dorothea-9606taxId_GeneCodis4.tsv")
TFregulonsNumbers <- table(TFregulons$synonyms)

resFiles <- list.files("random_list_analysis",pattern = "*.tsv",full.names = T)
resFile <- resFiles[1]
corrStats <- c("annotDB","nTFsRand","hypCorrFreqMeanTFs","walCorrFreqMeanTFs")
for (resFile in resFiles){
  resInfo <- unlist(strsplit(basename(resFile),split = "_"))
  annotDB <- resInfo[1]
  annotFile <- list.files(".",pattern = tolower(annotDB),full.names = T)
  nTFsRand <- gsub(".tsv","",resInfo[length(resInfo)])
  annots <- read.delim(annotFile)
  res <- read.delim(resFile)
  resHyp <- res[res$typeOfAnalysis == "Target_Hypergeom",]
  resWal <- res[res$typeOfAnalysis == "Target_NonCentral",]
  
  resHyp <- resHyp[order(resHyp$annotation_id),]
  resWal <- resWal[order(resWal$annotation_id),]
  
  annots <- annots[annots$annotation_id %in% resHyp$annotation_id,]
  
  meanNtfsDist <- c()
  for (annotid in resHyp$annotation_id){
    annotGenes <- annots$synonyms[annots$annotation_id == annotid]
    meanNtfs <- mean(na.omit(as.numeric(TFregulonsNumbers[annotGenes])))
    meanNtfsDist <- c(meanNtfsDist,meanNtfs)
  }
  
  resHyp$meanNtfs <- meanNtfsDist
  resWal$meanNtfs <- meanNtfsDist
  resHyp <- resHyp[!is.na(resHyp$meanNtfs),]
  resWal <- resWal[!is.na(resWal$meanNtfs),]
  
  write.table(cbind(resHyp$annotation_id,resHyp$meanNtfs),paste0(annotDB,nTFsRand,"annotMeanTFsDist.tsv"),sep = "\t")

  hypCorrFreqMeanTFs <- cor(resHyp$meanNtfs,resHyp$p_times)
  walCorrFreqMeanTFs <- cor(resWal$meanNtfs,resWal$p_times)  
  corrStats <- rbind(corrStats,c(annotDB,nTFsRand,hypCorrFreqMeanTFs,walCorrFreqMeanTFs))
}

corrStats <- as.data.frame(corrStats)
write.table(corrStats,"corrStats.tsv",sep = "\t")


