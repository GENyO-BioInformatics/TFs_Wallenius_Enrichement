## Install Packages

requireNamespace("BiocManager", quietly = TRUE) 
install.packages("BiocManager")
BiocManager::install('OmnipathR')

library(OmnipathR)
library(ggVennDiagram)
library(ggplot2)

##############################
##############################

## Fetch Data Clean and Save TF-Target Data
# dorotheaTFsGRN <- OmnipathR::dorothea(organism=9606, genesymbols=TRUE, loops=TRUE)
collectriTFsGRN <- OmnipathR::collectri(organism=9606, genesymbols=TRUE, loops=TRUE)
collectriTFsGRN <- as.data.frame(collectriTFsGRN)

collectriTFsGRN$source_genesymbol ## TF Gene
collectriTFsGRN$target_genesymbol ## Target Gene 
collectriTFsGRN_clean <- data.frame(tf=collectriTFsGRN$source_genesymbol,target=collectriTFsGRN$target_genesymbol,confidence='A',org='9606')

write.table(collectriTFsGRN, file = "data/collectri_raw.tsv",sep = '\t',row.names = F,col.names = T)
write.table(collectriTFsGRN_clean, file = "data/collectri.tsv",sep = '\t',row.names = F,col.names = T)

##############################
##############################

## Fetch GO Annotation Data 

# OMNIPATH
GObp <- OmnipathR::go_annot_download('human',aspects = 'P')
# GeneCodis4
gcGObp <- read.delim('data/GO_BP.tsv',header = T,sep = '\t')
gcGObp <- gcGObp[gcGObp$organism == 9606,]

## COMPARE GO data OMNIPATH and GeneCodis4
omnipathGOgenePairs <- paste(GObp$db_object_symbol,GObp$go_id,sep = '-')
gc4GOgenePairs <-paste(gcGObp$symbol,gcGObp$annotation_id,sep = '-')

nameOut <- "data/GOgenePairs.png"
toVennList <- list(omnipathGOgenePairs = omnipathGOgenePairs, gc4GOgenePairs = gc4GOgenePairs)
ggvenn <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
          scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
          theme(legend.position = "none")
ggsave(filename = nameOut, plot = ggvenn, width = 30, height = 15, units = 'cm', dpi = 'print')

nameOut <- "data/GOs.png"
toVennList <- list(omnipathGOs = GObp$go_id, gc4GOs = gcGObp$annotation_id)
ggvenn <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
ggsave(filename = nameOut, plot = ggvenn, width = 30, height = 15, units = 'cm', dpi = 'print')

nameOut <- "data/genesOfGO.png"
toVennList <- list(omnipathGenesOfGOs = GObp$db_object_symbol, gc4GenesOfGOs = gcGObp$symbol)
ggvenn <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
ggsave(filename = nameOut, plot = ggvenn, width = 30, height = 15, units = 'cm', dpi = 'print')


## Fetch KEGG Annotation Data 

KEGG <- OmnipathR::kegg_pathway_annotations()
gcKEGG <- read.delim('data/KEGG.tsv',header = T,sep = '\t')
gcKEGG <- gcKEGG[gcKEGG$organism == 9606,]

## COMPARE KEGG data OMNIPATH and GeneCodis4
omnipathKEGGgenePairs <- paste(KEGG$genesymbol,KEGG$pathway_id,sep = '-')
gc4KEGGgenePairs <-paste(gcKEGG$symbol,gcKEGG$annotation_id,sep = '-')

nameOut <- "data/KEGGgenePairs.png"
toVennList <- list(omnipathKEGGgenePairs = omnipathKEGGgenePairs, gc4KEGGgenePairs = gc4KEGGgenePairs)
ggvenn <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
ggsave(filename = nameOut, plot = ggvenn, width = 30, height = 15, units = 'cm', dpi = 'print')

nameOut <- "data/KEGGids.png"
toVennList <- list(omnipathKEGGs = KEGG$pathway_id, gc4KEGGs = gcKEGG$annotation_id)
ggvenn <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
ggsave(filename = nameOut, plot = ggvenn, width = 30, height = 15, units = 'cm', dpi = 'print')

nameOut <- "data/genesOfKEGG.png"
toVennList <- list(omnipathGenesOfKEGGs = KEGG$genesymbol, gc4GenesOfKEGGs = gcKEGG$symbol)
ggvenn <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
ggsave(filename = nameOut, plot = ggvenn, width = 30, height = 15, units = 'cm', dpi = 'print')

##############################
##############################
## CollecTri Annotations Gene Convergence Study

collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')

gcGObp <- read.delim('data/GO_BP.tsv',header = T,sep = '\t')
gcGObp <- gcGObp[gcGObp$organism == 9606,]

gcKEGG <- read.delim('data/KEGG.tsv',header = T,sep = '\t')
gcKEGG <- gcKEGG[gcKEGG$organism == 9606,]

gcReactome <- read.delim('data/Reactome.tsv',header = T,sep = '\t')
gcReactome <- gcReactome[gcReactome$organism == 9606,]

gcWikiPathways <- read.delim('data/WikiPathways.tsv',header = T,sep = '\t')
gcWikiPathways <- gcWikiPathways[gcWikiPathways$organism == 9606,]

nameOut <- "data/symbolsCollectriTFsnAnnotations.png"
toVennList <- list(GO=gcGObp$symbol, KEGG=gcKEGG$symbol, 
                   Reactome=gcReactome$symbol, WikiPathways=gcWikiPathways$symbol,
                   collectri=collectriTFsGRN$source_genesymbol)
ggvenn <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
ggsave(filename = nameOut, plot = ggvenn, width = 30, height = 15, units = 'cm', dpi = 'print')

nameOut <- "data/symbolsCollectriTargetsnAnnotations.png"
toVennList <- list(GO=gcGObp$symbol, KEGG=gcKEGG$symbol, 
                   Reactome=gcReactome$symbol, WikiPathways=gcWikiPathways$symbol,
                   collectri=collectriTFsGRN$target_genesymbol)
ggvenn <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
ggsave(filename = nameOut, plot = ggvenn, width = 30, height = 15, units = 'cm', dpi = 'print')

## CollecTri Curation Effort Average By Annotation
collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv'); annotationFile <- annotationFiles[1]
annotationsCurationEffortStats <- c()
for (annotationFile in annotationFiles){
  annotation <- gsub('.tsv','',basename(annotationFile))
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == 9606,]
  
  collectri_annTargets <- collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% annotationDF$symbol,]
  annotationDF_Targets <- annotationDF[annotationDF$symbol %in% collectri_annTargets$target_genesymbol,]
  
  collectriAnnotMerge <- merge(collectri_annTargets, annotationDF_Targets, by.x="target_genesymbol", by.y = "symbol")
  collectriAnnotMerge <- collectriAnnotMerge[,c('annotation_id','target_genesymbol','source_genesymbol','curation_effort')]
  collectriAnnotMerge <- collectriAnnotMerge[!duplicated(collectriAnnotMerge),]
  
  collectriAnnotMeanCEffort <- aggregate(collectriAnnotMerge$curation_effort, by = list(collectriAnnotMerge$annotation_id), FUN = mean)
  collectriAnnotMedianCEffort <- aggregate(collectriAnnotMerge$curation_effort, by = list(collectriAnnotMerge$annotation_id), FUN = median)
  collectriAnnotSumCEffort <- aggregate(collectriAnnotMerge$curation_effort, by = list(collectriAnnotMerge$annotation_id), FUN = sum)
  colnames(collectriAnnotMeanCEffort) <- c("annotation_id", "targets_average_curation_effort")
  colnames(collectriAnnotMedianCEffort) <- c("annotation_id", "targets_median_curation_effort")
  colnames(collectriAnnotSumCEffort) <- c("annotation_id", "targets_sum_curation_effort")
  collectriAnnotCEffortStats <- merge(merge(collectriAnnotMeanCEffort,collectriAnnotMedianCEffort),collectriAnnotSumCEffort)

  collectri_annTFs <- collectriTFsGRN[collectriTFsGRN$source_genesymbol %in% annotationDF$symbol,]
  annotationDF_TFs <- annotationDF[annotationDF$symbol %in% collectri_annTFs$target_genesymbol,]
  
  TFs_collectriAnnotMerge <- merge(collectri_annTFs, annotationDF_TFs, by.x="source_genesymbol", by.y = "symbol")
  TFs_collectriAnnotMerge <- TFs_collectriAnnotMerge[,c('annotation_id','target_genesymbol','source_genesymbol','curation_effort')]
  TFs_collectriAnnotMerge <- TFs_collectriAnnotMerge[!duplicated(TFs_collectriAnnotMerge),]
  
  TFs_collectriAnnotMeanCEffort <- aggregate(TFs_collectriAnnotMerge$curation_effort, by = list(TFs_collectriAnnotMerge$annotation_id), FUN = mean)
  TFs_collectriAnnotMedianCEffort <- aggregate(TFs_collectriAnnotMerge$curation_effort, by = list(TFs_collectriAnnotMerge$annotation_id), FUN = median)
  TFs_collectriAnnotSumCEffort <- aggregate(TFs_collectriAnnotMerge$curation_effort, by = list(TFs_collectriAnnotMerge$annotation_id), FUN = sum)
  colnames(TFs_collectriAnnotMeanCEffort) <- c("annotation_id", "tfs_average_curation_effort")
  colnames(TFs_collectriAnnotMedianCEffort) <- c("annotation_id", "tfs_median_curation_effort")
  colnames(TFs_collectriAnnotSumCEffort) <- c("annotation_id", "tfs_sum_curation_effort")
  TFs_collectriAnnotCEffortStats <- merge(merge(TFs_collectriAnnotMeanCEffort,TFs_collectriAnnotMedianCEffort),TFs_collectriAnnotSumCEffort)
  
  collectriAnnotCEffortStatsFull <-  merge(collectriAnnotCEffortStats, TFs_collectriAnnotCEffortStats,all = T)  
  collectriAnnotCEffortStatsFull$annotation <- annotation
  annotationsCurationEffortStats <- rbind(annotationsCurationEffortStats,collectriAnnotCEffortStatsFull)
}
write.table(annotationsCurationEffortStats, file = 'data/observations/annotationsCurationEffortStats.tsv' ,sep = '\t',
            row.names = F, col.names = T)

## Top Curated Effort Annotations
collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv'); annotationFile <- annotationFiles[1]

topKnown <- 50
for (annotationFile in annotationFiles){
  annotation <- gsub('.tsv','',basename(annotationFile))
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == 9606,]
  cat("Annotation:", annotation, "\n")
  for (top in (topKnown-20):topKnown){
    curationEffortCutOff <- unique(collectriTFsGRN$curation_effort[order(collectriTFsGRN$curation_effort,decreasing = T)])[top]
    mostKnownTargets <- unique(collectriTFsGRN$target_genesymbol[collectriTFsGRN$curation_effort >= curationEffortCutOff])
    mostKnownTFs <- unique(collectriTFsGRN$source_genesymbol[collectriTFsGRN$curation_effort >= curationEffortCutOff])
    annotationsCurationEffortStats <- c()
    # cat("Top:", top,
    #     "\tTarget annotations", 
    #     length(unique(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTargets])),
    #     "\tand\tTFs annotations", 
    #     length(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTFs]),
    #     "\n"
    # )
    cat("Top: ", top,
        "\tTarget most common annotation: ",
        max(table(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTargets])),
        "/",
        mean(table(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTargets])),
        "\tTFs most common annotation: ",
        max(table(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTFs])),
        "/",
        mean(table(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTFs])),
        "\n")
  }
}

## Mann-Whitney U test (Wilcoxon rank-sum test)
annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv'); annotationFile <- annotationFiles[1]
resultsDF <- c()
for (annotationFile in annotationFiles){
  annotation <- gsub('.tsv','',basename(annotationFile))
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == 9606,]
  collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
  collectriTFsGRN <- collectriTFsGRN[!grepl("_",collectriTFsGRN$source),]
  
  collectriTFsGRN <- collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% annotationDF$symbol, ]
  annotationDF <- annotationDF[annotationDF$symbol %in% collectriTFsGRN$target_genesymbol,]
  
  mergedData <- merge(annotationDF,collectriTFsGRN,by.x = "symbol", by.y="target_genesymbol")
  
  TFsperTargetAnnot <- as.data.frame(t(table(collectriTFsGRN$target_genesymbol)))[,2:3]
  colnames(TFsperTargetAnnot)[1] <- "target_genesymbol"
  
  for (annot in unique(annotationDF$annotation_id)){
    # annot <- unique(annotationDF$annotation_id)[1]
    TargetsAnnot <- unique(mergedData$symbol[mergedData$annotation_id == annot])
    annotTargetsPerTF <- TFsperTargetAnnot$Freq[TFsperTargetAnnot$target_genesymbol %in% TargetsAnnot]
    res <- wilcox.test(annotTargetsPerTF, TFsperTargetAnnot$Freq)
    resultsDF <- rbind(resultsDF,cbind(annot,res$p.value,annotation))
  }
}  
resultsDF <- as.data.frame(resultsDF)
colnames(resultsDF) <- c("annotation_id","wilcox.test","annotation")
write.table(resultsDF, file = 'data/observations/wilcox.test.tsv' ,sep = '\t',
            row.names = F, col.names = T)


## CollecTri Annotations Coverage By Evidence Study

#getAnnotationCoverage <- function(annotationFile,evidence,regulon=collectriTFsGRN){
evidences <- c('n_references','curation_effort'); evidence <- evidences[1]
annotationsCoveragePerEvidence <- c()
for (evidence in evidences){
  for (annotationFile in annotationFiles){
    annotation <- gsub('.tsv','',basename(annotationFile))
    outfile <- gsub('.tsv',paste0('_CollectriCoverage_',evidence,'.tsv'),annotationFile)
    outfile <- gsub('data','data/observations',outfile)
  
    annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
    annotationDF <- annotationDF[annotationDF$organism == 9606,]
    
    totalAnnotations <- length(unique(annotationDF$annotation_id))
    for (evidenceScore in sort(unique(collectriTFsGRN[,evidence]))){
      selection <- which(collectriTFsGRN[,evidence] >= evidenceScore)
      TFsAccepted <- unique(collectriTFsGRN$source_genesymbol[selection])
      TargetsAccepted <- unique(collectriTFsGRN$target_genesymbol[selection])
      nTFsAccepted <- length(TFsAccepted)
      nTargetsAccepted <- length(TargetsAccepted)
      TFsAnnotationsInEvidence <- length(unique(annotationDF$annotation_id[annotationDF$symbol %in% TFsAccepted]))
      TargetsAnnotationsInEvidence <- length(unique(annotationDF$annotation_id[annotationDF$symbol %in% TargetsAccepted]))
      propor_TFsAnnotationsInEvidence <- TFsAnnotationsInEvidence / totalAnnotations * 100
      propot_TargetsAnnotationsInEvidence <- TargetsAnnotationsInEvidence / totalAnnotations * 100
      
      annotationsCoveragePerEvidence <- rbind(annotationsCoveragePerEvidence,
                                             cbind(annotation=annotation,
                                                   evidence = evidence,
                                                   minEvidenceScore = evidenceScore,
                                                   nTFsAccepted = nTFsAccepted, 
                                                   TFsAnnotationsInEvidence = TFsAnnotationsInEvidence,
                                                   propor_TFsAnnotationsInEvidence = propor_TFsAnnotationsInEvidence,
                                                   nTargetsAccepted = nTargetsAccepted,
                                                   TargetsAnnotationsInEvidence = TargetsAnnotationsInEvidence,
                                                   propor_TargetsAnnotationsInEvidence= propot_TargetsAnnotationsInEvidence))
    }
  }
}
write.table(annotationsCoveragePerEvidence, file = 'data/observations/annotationsCoveragePerEvidence.tsv' ,sep = '\t',
            row.names = F, col.names = T)
View(annotationsCoveragePerEvidence)

##############################
##############################
## CollecTri TF - Target Distribution - Per Annotation DB

collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv'); annotationFile <- annotationFiles[1]
TFsPerTargetdf_full <- c()
for (annotationFile in annotationFiles){
  annotation <- gsub('.tsv','',basename(annotationFile))
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == 9606,]
  collectriTFsGRN_filtered <- collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% annotationDF$symbol,]
  TFsPerTargetdf <- as.data.frame(table(collectriTFsGRN_filtered$target_genesymbol))
  colnames(TFsPerTargetdf) <- c('target','numberOfTFsAssociated')
  TFsPerTargetdf$db <- annotation
  TFsPerTargetdf_full <- rbind(TFsPerTargetdf_full, TFsPerTargetdf)
}
View(TFsPerTargetdf_full)

write.table(TFsPerTargetdf_full, file = 'data/observations/TFsPerTargetdf.tsv',
            sep = '\t', row.names = F, col.names = T)


quantile(TFsPerTargetdf$numberOfTFsAssociated)
mean(TFsPerTargetdf$numberOfTFsAssociated)
median(TFsPerTargetdf$numberOfTFsAssociated)
ux <- unique(TFsPerTargetdf$numberOfTFsAssociated)
tab <- tabulate(match(TFsPerTargetdf$numberOfTFsAssociated, ux))
ux[tab == max(tab)]

##############################
##############################
# WHAT ABOUT THE TF COMPLEXES? Shall we remove them?
collectriTFsGRN[!grepl('COMPLEX',collectriTFsGRN$source),]
sum(grepl('COMPLEX',collectriTFsGRN$source))
# NO because we retrieve the targets where there are no COMPLEXES
sum(grepl('COMPLEX',collectriTFsGRN$target))


##############################
##############################
## Generate random TFs lists
sizes <- c(3, 10, 15, 20); size <- sizes[1]

set.seed(9)
for (size in sizes){
  for (n in 1:1000){
    outdir <- paste0('data/TFsLists/size',size)
    dir.create(outdir,recursive = T,showWarnings = F)
    TFsList <- sample(collectriTFsGRN$source_genesymbol,size=size,replace = F)
    write.table(TFsList, file = paste0(outdir,'/',n,'.txt'),
                sep = '\t', row.names = F, col.names = F,quote = F)
  }
}

##############################
##############################
## Launch
'cmder_analyse_random_TFlists.py'

#############################
#############################

# Get The Average of TFs per Annotations
collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
source2TargetDF <- collectriTFsGRN[,c('source_genesymbol','target_genesymbol')]
collectriTFsGRN <- collectriTFsGRN[!duplicated(source2TargetDF),]
View(collectriTFsGRN[duplicated(source2TargetDF),])

annotationDBs <- c('GO_BP','KEGG','Reactome','WikiPathways')
statsDF <- c()
# annotationDB <- 'Reactome'; annotation <- 'R-HSA-9646303' 
for (annotationDB in annotationDBs){
  annotationFile <- paste0('data/',annotationDB,'.tsv')
  annotationDBdf <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDBdf <- annotationDBdf[annotationDBdf$organism == 9606,]
  annotations <- unique(annotationDBdf$annotation_id)
  for (annotation in annotations){
    annotationGenes <- annotationDBdf$symbol[annotationDBdf$annotation_id == annotation] 
    targetsInAnnotation <- unique(collectriTFsGRN$target_genesymbol[collectriTFsGRN$target_genesymbol %in% annotationGenes])
    TFsInAnnotation <- unique(collectriTFsGRN$source_genesymbol[collectriTFsGRN$source_genesymbol %in% annotationGenes])
    tfsAssociatedWithTargetsAnnotation <- unique(collectriTFsGRN$source_genesymbol[collectriTFsGRN$target_genesymbol %in% targetsInAnnotation])
    averageTFsPerTargetInAnnotation <- mean(table(collectriTFsGRN$target_genesymbol[collectriTFsGRN$target_genesymbol %in% targetsInAnnotation]))
    averageTargetsPerTFInAnnotation <- mean(table(collectriTFsGRN$source_genesymbol[collectriTFsGRN$target_genesymbol %in% targetsInAnnotation]))
    statsPerAnnotation <- c(annotationDB,annotation, length(annotationGenes), length(TFsInAnnotation), length(targetsInAnnotation), length(tfsAssociatedWithTargetsAnnotation), averageTFsPerTargetInAnnotation, averageTargetsPerTFInAnnotation)
    statsDF <- rbind(statsDF,statsPerAnnotation)
  }
}

View(collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% targetsInAnnotation, c('source_genesymbol','target_genesymbol')])

# 'annotationDB' -> The annotation database
# 'annotation' -> The annotation
# 'annotationGenes' -> The genes associated to annotation in Functional DB
# 'TFsInAnnotation' -> The number of annotationGenes that are TFs
# 'targetsInAnnotation' ->  The annotationGenes that are targets according to Collectri 
# 'tfsAssociatedWithTargetsAnnotation' -> The TFs that point to targetsInAnnotation
# 'averageTFsPerTargetInAnnotation' -> The mean number of TFs per targetsInAnnotation
# 'averageTargetsPerTFInAnnotation' -> The mean number of targets per TF associated (tfsAssociatedWithTargetsAnnotation)

statsDF <- as.data.frame(statsDF)
colNames <- c('annotationDB','annotation','annotationGenes','TFsInAnnotation','targetsInAnnotation','tfsAssociatedWithTargetsAnnotation','averageTFsPerTargetInAnnotation','averageTargetsPerTFInAnnotation')
colnames(statsDF) <- colNames
statsDF <- type.convert(statsDF,as.is=T)
statsDF$averageTFsPerTargetInAnnotation[is.nan(statsDF$averageTFsPerTargetInAnnotation)] <- 0
statsDF$averageTargetsPerTFInAnnotation[is.nan(statsDF$averageTargetsPerTFInAnnotation)] <- 0
View(statsDF)

write.table(statsDF, file = 'data/observations/TFsNtargetsPerAnnotationDistribution.tsv',
            sep = '\t', row.names = F, col.names = T, quote = F)


averageTFsPerAnnotation
#table(collectriTFsGRN$source_genesymbol[collectriTFsGRN$target_genesymbol %in% targetsInAnnotation])
#collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% targetsInAnnotation,c('source_genesymbol','target_genesymbol')]

collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')

topCuratedScore <- unique(collectriTFsGRN$curation_effort)[order(unique(collectriTFsGRN$curation_effort),decreasing = T)[5]]
topCuratedScore <- 200
collectriTFsGRN_top5Curated <- collectriTFsGRN[collectriTFsGRN$curation_effort >= topCuratedScore,]

library(igraph)

toIgraph <- collectriTFsGRN_top5Curated[c("source_genesymbol","target_genesymbol","is_stimulation","curation_effort")]
igraphObj <- graph_from_data_frame(toIgraph[!grepl("_",toIgraph$source_genesymbol),], directed = T, vertices = NULL)

png("data/observations/networkAbove200curationeffort.png",width = 10, height = 10, units = "cm", res = 300)
plot(igraphObj,  edge.arrow.size=.4, edge.color="grey",
     vertex.color="orange", vertex.frame.color="#ffffff",
     vertex.label=V(igraphObj)$media, vertex.label.color="black",
     vertex.shape="none",layout=layout_with_kk)
dev.off()
