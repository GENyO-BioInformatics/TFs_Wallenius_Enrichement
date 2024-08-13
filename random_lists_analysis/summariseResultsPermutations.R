library(vroom)
library(data.table)
library(tidyverse)
npermutations <- 1000

# PID 3486739 <- 

resultsFiles <- Sys.glob('data/EnrResults/*/*/*[0-9].tsv', dirmark = FALSE)
annotNsize <- do.call(rbind, strsplit(resultsFiles,split = '/'))[,3:4]
annotNsize <- annotNsize[!duplicated(annotNsize),]
idx <- 1
for (idx in 1:nrow(annotNsize)){
  db <- annotNsize[idx,1]
  size <- annotNsize[idx,2]
  resultsFiles <- Sys.glob(file.path('data/EnrResults',db,size,'*[0-9].tsv'), dirmark = FALSE)
  resultsFile <- resultsFiles[1]; resultsFile; resultsFiles[2]
  myColumns <- c('terms',
                 'TFs_Hypergeom','ranking_TFs_Hypergeom',
                 'Target_Hypergeom','ranking_Target_Hypergeom',
                 'Target_NonCentral','ranking_Target_NonCentral', 
                 'file')
  # ÑAPA 1
  myColumns <- c(myColumns, c('TFsTargetsUni_Hypergeom','ranking_TFsTargetsUni_Hypergeom')) 
  # FIN ÑAPA 1
  
  masterDF <- setNames(data.frame(matrix(ncol = length(myColumns), nrow = 0)), myColumns)
  for (resultsFile in resultsFiles){
    #resultsDF <- vroom(resultsFile,show_col_types = FALSE)
    # ÑAPA 2
    oldresultsDF <- vroom(resultsFile,show_col_types = FALSE)
    newResFile <- gsub(".tsv","_TFsTargetUniv.tsv",resultsFile)
    newresultsDF <- vroom(newResFile,show_col_types = FALSE)
    colnames(newresultsDF) <- gsub("TFs","TFsTargetsUni",colnames(newresultsDF))
    resultsDF <- merge(oldresultsDF,newresultsDF)
    # FIN ÑAPA 2
    
    resultsDF[order(resultsDF$TFsTargetsUni_Hypergeom,resultsDF$terms),'ranking_TFsTargetsUni_Hypergeom'] <- c(1:nrow(resultsDF))
    resultsDF[order(resultsDF$TFs_Hypergeom,resultsDF$terms),'ranking_TFs_Hypergeom'] <- c(1:nrow(resultsDF))
    resultsDF[order(resultsDF$Target_Hypergeom,resultsDF$terms),'ranking_Target_Hypergeom'] <- c(1:nrow(resultsDF))
    resultsDF[order(resultsDF$Target_NonCentral,resultsDF$terms),'ranking_Target_NonCentral'] <- c(1:nrow(resultsDF))
    resultsDF$file <- basename(resultsFile)
    masterDF <- rbind(masterDF,resultsDF[,myColumns])
  }
  allResultsOutFile <- paste0(dirname(resultsFile[1]),'.tsv')
  write.table(masterDF, allResultsOutFile, sep='\t',quote = F,row.names = F,col.names = T)
  
  TFs_HypergeomCount <- aggregate(TFs_Hypergeom ~ terms, masterDF, function(x) sum(x <= 0.05) / npermutations)
  colnames(TFs_HypergeomCount)[2] <- 'TFs_HypergeomCount'
  Target_HypergeomCount <- aggregate(Target_Hypergeom ~ terms, masterDF, function(x) sum(x <= 0.05) / npermutations)
  colnames(Target_HypergeomCount)[2] <- 'Target_HypergeomCount'
  Target_NonCentralCount <- aggregate(Target_NonCentral ~ terms, masterDF, function(x) sum(x <= 0.05) / npermutations)
  colnames(Target_NonCentralCount)[2] <- 'Target_NonCentralCount'
  
  # ÑAPA 3.1
  TFsTargetsUni_HypergeomCount <- aggregate(TFsTargetsUni_Hypergeom ~ terms, masterDF, function(x) sum(x <= 0.05) / npermutations)
  colnames(TFsTargetsUni_HypergeomCount)[2] <- 'TFsTargetsUni_HypergeomCount'
  ######
  
  ranking_TFs_Hypergeom_Median <- aggregate(ranking_TFs_Hypergeom ~ terms, masterDF, function(x) median(x)) 
  colnames(ranking_TFs_Hypergeom_Median)[2] <- 'ranking_TFs_Hypergeom_Median'
  ranking_Target_Hypergeom_Median <- aggregate(ranking_Target_Hypergeom ~ terms, masterDF, function(x) median(x)) 
  colnames(ranking_Target_Hypergeom_Median)[2] <- 'ranking_Target_Hypergeom_Median'
  ranking_Target_NonCentral_Median <- aggregate(ranking_Target_NonCentral ~ terms, masterDF, function(x) median(x)) 
  colnames(ranking_Target_NonCentral_Median)[2] <- 'ranking_Target_NonCentral_Median'
  
  # ÑAPA 3.2
  ranking_TFsTargetsUni_Hypergeom_Median <- aggregate(ranking_TFsTargetsUni_Hypergeom ~ terms, masterDF, function(x) median(x)) 
  colnames(ranking_TFsTargetsUni_Hypergeom_Median)[2] <- 'ranking_TFsTargetsUni_Hypergeom_Median'
  ######
  
  # ranking_TFs_Hypergeom_IQR <- aggregate(ranking_TFs_Hypergeom ~ terms, masterDF, function(x) IQR(x)) 
  # colnames(ranking_TFs_Hypergeom_IQR)[2] <- 'ranking_TFs_Hypergeom_IQR'
  # ranking_Target_Hypergeom_IQR <- aggregate(ranking_Target_Hypergeom ~ terms, masterDF, function(x) IQR(x)) 
  # colnames(ranking_Target_Hypergeom_IQR)[2] <- 'ranking_Target_Hypergeom_IQR'
  # ranking_Target_NonCentral_IQR <- aggregate(ranking_Target_NonCentral ~ terms, masterDF, function(x) IQR(x)) 
  # colnames(ranking_Target_NonCentral_IQR)[2] <- 'ranking_Target_NonCentral_IQR'

  ranking_TFs_Hypergeom_Q1 <- aggregate(ranking_TFs_Hypergeom ~ terms, masterDF, function(x) quantile(x, prob=.25, type=1))
  colnames(ranking_TFs_Hypergeom_Q1)[2] <- 'ranking_TFs_Hypergeom_Q1'
  ranking_Target_Hypergeom_Q1 <- aggregate(ranking_Target_Hypergeom ~ terms, masterDF, function(x) quantile(x, prob=.25, type=1))
  colnames(ranking_Target_Hypergeom_Q1)[2] <- 'ranking_Target_Hypergeom_Q1'
  ranking_Target_NonCentral_Q1 <- aggregate(ranking_Target_NonCentral ~ terms, masterDF, function(x) quantile(x, prob=.25, type=1))
  colnames(ranking_Target_NonCentral_Q1)[2] <- 'ranking_Target_NonCentral_Q1'

  # ÑAPA 3.3
  ranking_TFsTargetsUni_Hypergeom_Q1 <- aggregate(ranking_TFsTargetsUni_Hypergeom ~ terms, masterDF, function(x) quantile(x, prob=.25, type=1)) 
  colnames(ranking_TFsTargetsUni_Hypergeom_Q1)[2] <- 'ranking_TFsTargetsUni_Hypergeom_Q1'
  #######

  ranking_TFs_Hypergeom_Q3 <- aggregate(ranking_TFs_Hypergeom ~ terms, masterDF, function(x) quantile(x, prob=.75, type=1))
  colnames(ranking_TFs_Hypergeom_Q3)[2] <- 'ranking_TFs_Hypergeom_Q3'
  ranking_Target_Hypergeom_Q3 <- aggregate(ranking_Target_Hypergeom ~ terms, masterDF, function(x) quantile(x, prob=.75, type=1))
  colnames(ranking_Target_Hypergeom_Q3)[2] <- 'ranking_Target_Hypergeom_Q3'
  ranking_Target_NonCentral_Q3 <- aggregate(ranking_Target_NonCentral ~ terms, masterDF, function(x) quantile(x, prob=.75, type=1))
  colnames(ranking_Target_NonCentral_Q3)[2] <- 'ranking_Target_NonCentral_Q3'

  # ÑAPA 3.4
  ranking_TFsTargetsUni_Hypergeom_Q3 <- aggregate(ranking_TFsTargetsUni_Hypergeom ~ terms, masterDF, function(x) quantile(x, prob=.75, type=1)) 
  colnames(ranking_TFsTargetsUni_Hypergeom_Q3)[2] <- 'ranking_TFsTargetsUni_Hypergeom_Q3'
  #######
  
  statisticsDF <- reduce(list(TFs_HypergeomCount,
                              ranking_TFs_Hypergeom_Median,
                              ranking_TFs_Hypergeom_Q1,
                              ranking_TFs_Hypergeom_Q3,
                              TFsTargetsUni_HypergeomCount,
                              ranking_TFsTargetsUni_Hypergeom_Median,
                              ranking_TFsTargetsUni_Hypergeom_Q1,
                              ranking_TFsTargetsUni_Hypergeom_Q3,
                              Target_HypergeomCount,
                              ranking_Target_Hypergeom_Median,
                              ranking_Target_Hypergeom_Q1,
                              ranking_Target_Hypergeom_Q3,
                              Target_NonCentralCount,
                              ranking_Target_NonCentral_Median,
                              ranking_Target_NonCentral_Q1,
                              ranking_Target_NonCentral_Q3),
                         full_join, by='terms')
  
  statisticsDF$size <- gsub('size','',size)
  statisticsDF$annotation <- db
  statisticsOutFile <- paste0(dirname(resultsFile[1]),'_statistics.tsv')
  write.table(statisticsDF, statisticsOutFile, sep='\t',quote = F,row.names = F,col.names = T)
}

summarystatsFiles <- list.files('data/EnrResults',pattern = 'statistics.tsv$', recursive = T,full.names = T)
masterSummaryDF <- read.delim(summarystatsFiles[1])
for (summarystatsFile in summarystatsFiles[-1]){
  masterSummaryDF <- rbind(masterSummaryDF,read.delim(summarystatsFile))
}
masterSummaryOutFile <- 'data/EnrResults/masterSummary.tsv'
write.table(masterSummaryDF, masterSummaryOutFile, sep='\t',quote = F,row.names = F,col.names = T)

# CREATION OF MELTED TABLES TO EASILY EXPLORE IN EXCEL #############
masterDF <- read.delim("data/EnrResults/masterSummary.tsv") 
masterDF$annotation <- gsub('_', ' ', masterDF$annotation)
annotInfo <- read.delim("data/annotation_info_table.tsv")
masterDF <- merge(annotInfo,masterDF,by.x ="annotation_id",by.y = "terms",all.y = T)

fullContent <- as.data.table(masterDF)

idCols <- c("annotation_id", "term","size","annotation")
permutationsCountCols <- colnames(fullContent)[grep("Count",colnames(fullContent))]
permutMeltedDF <- melt(fullContent, id.vars=idCols, measure.vars=permutationsCountCols)
colnames(permutMeltedDF) <- c(idCols,"typeOfAnalysis","p_times")
permutMeltedDF$typeOfAnalysis <- gsub("Count","",permutMeltedDF$typeOfAnalysis)

rankMedianCols <- colnames(fullContent)[grep("Median",colnames(fullContent))]
rankMedianMeltedDF <- melt(fullContent, id.vars=idCols, measure.vars=rankMedianCols)
colnames(rankMedianMeltedDF) <- c(idCols,"typeOfAnalysis","RankingMedian")
rankMedianMeltedDF$typeOfAnalysis <- gsub("ranking_|_Median","",rankMedianMeltedDF$typeOfAnalysis)

q1Cols <- colnames(fullContent)[grep("Q1",colnames(fullContent))]
q1MeltedDF <- melt(fullContent, id.vars=idCols, measure.vars=q1Cols)
colnames(q1MeltedDF) <- c(idCols,"typeOfAnalysis","RankingQ1")
q1MeltedDF$typeOfAnalysis <- gsub("ranking_|_Q1","",q1MeltedDF$typeOfAnalysis)

q3Cols <- colnames(fullContent)[grep("Q3",colnames(fullContent))]
q3MeltedDF <- melt(fullContent, id.vars=idCols, measure.vars=q3Cols)
colnames(q3MeltedDF) <- c(idCols,"typeOfAnalysis","RankingQ3")
q3MeltedDF$typeOfAnalysis <- gsub("ranking_|_Q3","",q3MeltedDF$typeOfAnalysis)

fullContentMelted <- merge(merge(merge(permutMeltedDF,rankMedianMeltedDF),q1MeltedDF),q3MeltedDF)

TFsSEAresults <- fullContentMelted[fullContentMelted$typeOfAnalysis == "TFs_Hypergeom" & !is.na(fullContentMelted$p_times),]
write.table(TFsSEAresults,"data/EnrResults/TFsSEAresultsMelted.tsv",sep = "\t",quote = F,row.names = F,col.names = T)

TargetsSEAresults <- fullContentMelted[fullContentMelted$typeOfAnalysis != "TFs_Hypergeom",]
write.table(TargetsSEAresults,"data/EnrResults/TargetsSEAresultsMelted.tsv",sep = "\t",quote = F,row.names = F,col.names = T)

