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
