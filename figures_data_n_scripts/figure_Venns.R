################ The following figures are VennDiagrams of the TFs/Targets genes
################ that are common bettween Collectri and the annotations dbs.
"
These plot answer the following question:
1- Which universe of TFs/Targets genes should be used to perform the null (montecarlo)
   simulations?
   We should use the intersection bewteen dbs to obtain the random lists and to perform the statistical tests
"

library(ggVennDiagram)
library(ggplot2)

setwd("/home/adriangarciamoreno/Desktop/TFs_Wallenius_Enrichement")
outdir <- "Figures/venns/"
setColors <- c("orange","red","blue","purple", "black")
names(setColors) <- c("Collectri","GOBP","KEGG","Reactome", "WikiPathways")


collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')

gcGObp <- read.delim('data/GO_BP.tsv',header = T,sep = '\t')
gcGObp <- gcGObp[gcGObp$organism == 9606,]

gcKEGG <- read.delim('data/KEGG.tsv',header = T,sep = '\t')
gcKEGG <- gcKEGG[gcKEGG$organism == 9606,]

gcReactome <- read.delim('data/Reactome.tsv',header = T,sep = '\t')
gcReactome <- gcReactome[gcReactome$organism == 9606,]

gcWikiPathways <- read.delim('data/WikiPathways.tsv',header = T,sep = '\t')
gcWikiPathways <- gcWikiPathways[gcWikiPathways$organism == 9606,]

####### VennDiagram of the TFs of Collectri and how these are annotated
symbolsCollectriTFsnAnnotations <- list(Collectri=collectriTFsGRN$source_genesymbol,
                                        GOBP=gcGObp$symbol, KEGG=gcKEGG$symbol, 
                                        Reactome=gcReactome$symbol, WikiPathways=gcWikiPathways$symbol)

nameOut <- file.path(outdir,"figure_TFs_in_DBs.png")
ggvenn <- ggVennDiagram(symbolsCollectriTFsnAnnotations, set_color = setColors, lwd = 0.7,
                        label="count",label_alpha=0) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")+
  labs(title = "All Transcription Factors in DBs")
ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')


####### VennDiagram of the Targets of Collectri and how these are annotated
symbolsCollectriTargetsnAnnotations <- list(Collectri=collectriTFsGRN$target_genesymbol,
                                            GOBP=gcGObp$symbol, KEGG=gcKEGG$symbol, 
                                            Reactome=gcReactome$symbol, WikiPathways=gcWikiPathways$symbol)
nameOut <- file.path(outdir,"figure_Targets_in_DBs.png")
ggvenn <- ggVennDiagram(symbolsCollectriTargetsnAnnotations, set_color = setColors, lwd = 0.7, 
                        label="count",label_alpha=0) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")+
  labs(title = "All Target Genes in DBs")
ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')


####### VennDiagrams of the Collectri TFs/Targets and how these are annotated separately
for (annotationDB in setdiff(names(symbolsCollectriTFsnAnnotations),'Collectri')){
  mysubset <- c("Collectri",annotationDB)
  subList <- symbolsCollectriTFsnAnnotations[mysubset]
  nameOut <- file.path(outdir,paste0("figure_TFs_Collectri_n_",annotationDB,".png"))
  ggvenn <- ggVennDiagram(subList, set_color = setColors[mysubset], lwd = 0.7, label="count",label_alpha=0) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none")+
    labs(title = paste0("All Transcription Factors in Collectri and",annotationDB))+
    scale_x_continuous(expand = expansion(mult = .2))
  ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')
  
  subList <- symbolsCollectriTargetsnAnnotations[mysubset]
  nameOut <- file.path(outdir,paste0("figure_Targets_Collectri_n_",annotationDB,".png"))
  ggvenn <- ggVennDiagram(subList, set_color = setColors[mysubset], lwd = 0.7, label="count",label_alpha=0) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none") +
    labs(title = paste0("All Transcription Factors Target Genes in Collectri and",annotationDB)) +
    scale_x_continuous(expand = expansion(mult = .2))
  ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')
}


####### VennDiagrams in the of Collectri Targets Universe representation in annotations
### We generate subsets of Collectri by only selecting the targets of Collectri that
### are annotated in each database. Then we plot the VennDiagrams of the resulting subsets
### at the TFs, Targets and Interaction levels.

collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
## ISSUE: some TFs of Collectri are also annotated in the functional dbs and we will be missing them
##        given that we are only considering the Collectri targets
all(collectriTFsGRN$source_genesymbol %in% collectriTFsGRN$target_genesymbol)
all(collectriTFsGRN$source_genesymbol %in% gcGObp$symbol)
all(collectriTFsGRN$source_genesymbol %in% gcKEGG$symbol)
all(collectriTFsGRN$source_genesymbol %in% gcReactome$symbol)
all(collectriTFsGRN$source_genesymbol %in% gcWikiPathways$symbol)
all(collectriTFsGRN$target_genesymbol %in% gcGObp$symbol)
all(collectriTFsGRN$target_genesymbol %in% gcKEGG$symbol)
all(collectriTFsGRN$target_genesymbol %in% gcReactome$symbol)
all(collectriTFsGRN$target_genesymbol %in% gcWikiPathways$symbol)


# CTFU == Collectri TFs Universe (to only compare the info of collectri based on the TFs that are annotated in DB)
CTFU_TFsRepresentation <- TFsRepresentation <- list()
CTFU_TargetsRepresentation <- TargetsRepresentation <- list()
CTFU_InteractionsRepresentation <- InteractionsRepresentation <- list()

CTFU_TFsRepresentation[['Collectri']] <- TFsRepresentation[['Collectri']] <- unique(collectriTFsGRN$source_genesymbol)
CTFU_TargetsRepresentation[['Collectri']] <- TargetsRepresentation[['Collectri']] <- unique(collectriTFsGRN$target_genesymbol)
CTFU_InteractionsRepresentation[['Collectri']] <- InteractionsRepresentation[['Collectri']] <- unique(paste0(collectriTFsGRN$source_genesymbol,collectriTFsGRN$target_genesymbol))
  

# https://stackoverflow.com/questions/68875752/how-to-edit-ggvenndiagram-intersection-fill-region

annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv')
annotationFile <- annotationFiles[1]
for (annotationFile in annotationFiles){
  annotation <- gsub('_','',gsub('.tsv','',basename(annotationFile)))
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == 9606,]
  
  subCollectri <- collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% annotationDF$symbol,]
  TFsRepresentation[[annotation]] <- unique(subCollectri$source_genesymbol)
  TargetsRepresentation[[annotation]] <- unique(subCollectri$target_genesymbol)
  InteractionsRepresentation[[annotation]] <- unique(paste0(subCollectri$source_genesymbol,subCollectri$target_genesymbol))
  
  subCollectri <- collectriTFsGRN[collectriTFsGRN$source_genesymbol %in% annotationDF$symbol,]
  CTFU_TFsRepresentation[[annotation]] <- unique(subCollectri$source_genesymbol)
  CTFU_TargetsRepresentation[[annotation]] <- unique(subCollectri$target_genesymbol)
  CTFU_InteractionsRepresentation[[annotation]] <- unique(paste0(subCollectri$source_genesymbol,subCollectri$target_genesymbol))

}

nameOut <- file.path(outdir,paste0("figure_Collectri_TFs_Universe_in_DBs.png"))
ggvenn <- ggVennDiagram(TFsRepresentation, set_color = setColors, lwd = 0.7, label="count",label_alpha=0) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Collectri Transcription Factors Universe in DBs") +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')

nameOut <- file.path(outdir,paste0("figure_Collectri_Targets_Universe_in_DBs.png"))
ggvenn <- ggVennDiagram(TargetsRepresentation, set_color = setColors, lwd = 0.7, label="count",label_alpha=0) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Collectri Target Genes Universe in DBs") +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')

nameOut <- file.path(outdir,paste0("figure_Collectri_Interactions_Universe_in_DBs.png"))
ggvenn <- ggVennDiagram(InteractionsRepresentation, set_color = setColors, lwd = 0.7, label="count",label_alpha=0) +   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Collectri TF-Target Genes Interactions Universe in DBs") +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')


# CTFU
nameOut <- file.path(outdir,paste0("figure_CTFU_TFs_Universe_in_DBs.png"))
ggvenn <- ggVennDiagram(CTFU_TFsRepresentation, set_color = setColors, lwd = 0.7, label="count",label_alpha=0) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "CTFU Transcription Factors Universe in DBs") +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')

nameOut <- file.path(outdir,paste0("figure_CTFU_Targets_Universe_in_DBs.png"))
ggvenn <- ggVennDiagram(CTFU_TargetsRepresentation, set_color = setColors, lwd = 0.7, label="count",label_alpha=0) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "CTFU Target Genes Universe in DBs") +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')

nameOut <- file.path(outdir,paste0("figure_CTFU_Interactions_Universe_in_DBs.png"))
ggvenn <- ggVennDiagram(CTFU_InteractionsRepresentation, set_color = setColors, lwd = 0.7, label="count",label_alpha=0) +   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "CTFU TF-Target Genes Interactions Universe in DBs") +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(filename = nameOut, plot = ggvenn, width = 15, height = 15, units = 'cm', dpi = 'print')



