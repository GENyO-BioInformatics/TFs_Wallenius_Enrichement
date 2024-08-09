############# Plots Collectri Targets and TFs distribution in Annotations DBs 
"
This plot answers the questions:
1- How is the distribution of targets per annotation? 
2- How is the distribution of TFs per annotation?
"

library(vroom)
library(ggplot2)
library(ggrepel)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library('ggpp')


TFsNtargets$db <- gsub('_',' ',TFsNtargets$db)
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")


TFsNtargetsAnnotaDist <- vroom("data/observations/TFsNtargetsPerAnnotationDistribution.tsv")

TFsNtargetsAnnotaDist$TFsInAnnotation_proportion <- TFsNtargetsAnnotaDist$TFsInAnnotation / TFsNtargetsAnnotaDist$annotationGenes 
TFsNtargetsAnnotaDist$targetsInAnnotation_proportion <- TFsNtargetsAnnotaDist$targetsInAnnotation / TFsNtargetsAnnotaDist$annotationGenes

TFsNtargetsAnnotaDist$annotationDB <- gsub('_',' ',TFsNtargetsAnnotaDist$annotationDB)

ggsave("Figures/TFsNtargetsPerAnnotationDistribution.tiff", plot = gg,units = "cm",height = 12, width = 20,dpi = 500)

# toCorrelationMatrix <- TFsNtargetsAnnotaDist[c(3:ncol(TFsNtargetsAnnotaDist))]
# View(cor(toCorrelationMatrix))
# View(toCorrelationMatrix)
# 
# library(GGally)
# ggpairs(toCorrelationMatrix, title="correlogram with ggpairs()",ggplot2::aes(colour=annotationDB)) 

TFsNtargetsAnnotaDist$TFsInAnnotation

