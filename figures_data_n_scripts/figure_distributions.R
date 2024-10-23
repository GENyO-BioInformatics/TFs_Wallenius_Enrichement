setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),".."))
getwd()
library(ggplot2)
library(gridExtra)
library(stringr)
library(smplot2)
library(ggrepel)
library(ggpubr)
library(vroom)
library(stringr)
library(tidyr)
swr = function(string, nwrap=50) {paste(strwrap(string, width=nwrap), collapse="\n")}; swr = Vectorize(swr)

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

dataFile <- 'data/observations/probabilities.csv'
annotationsSelectionProbs <- read.delim(dataFile,sep = ',')
annotationsSelectionProbs$probability
colnames(annotationsSelectionProbs)[1:2] <- c("annotation","annotation_id")

dataFile <- 'random_lists_analysis/masterSimsResults.tsv'
simulationsResultsDF <- read.delim(dataFile)
colnames(simulationsResultsDF)[1] <- "annotation_id"
simulationsResultsDF$annotation <- gsub("_"," ",simulationsResultsDF$annotation)

dataFile <- 'data/observations/annotationsCurationEffortStats.tsv'
annotationsCurationEffortStats <- read.delim(dataFile)
annotationsCurationEffortStats$annotation <- gsub("_"," ",annotationsCurationEffortStats$annotation)
nrow(annotationsCurationEffortStats)

dataFile <- 'data/observations/TFsNtargetsPerAnnotationDistribution.tsv'
TFsNtargetsPerAnnotDF <- read.delim(dataFile)
colnames(TFsNtargetsPerAnnotDF)[1:2] <- c("annotation","annotation_id")
TFsNtargetsPerAnnotDF$annotation <- gsub("_"," ",TFsNtargetsPerAnnotDF$annotation)

TFsNtargetsPerAnnotDF$tfsAssociatedWithTargetsAnnotation

table(annotationsSelectionProbs$annotation)
table(simulationsResultsDF$annotation)
table(annotationsCurationEffortStats$annotation)
table(TFsNtargetsPerAnnotDF$annotation)

fullStatsDF <- merge(annotationsCurationEffortStats,
                     merge(TFsNtargetsPerAnnotDF,
                      merge(annotationsSelectionProbs,simulationsResultsDF)))

annotationsDistributionsCol <- c("pTargetsF","pTargetsW", "probability")
# The following variables are intereseting: 
# "tfsAssociatedWithTargetsAnnotation" shows an almost identicial distribution as "probability" ; 
# "averageTargetsPerTFInAnnotation" distort the plot
# annotationsDistributions$averageTargetsPerTFInAnnotation <- annotationsDistributions$averageTargetsPerTFInAnnotation / max(annotationsDistributions$averageTargetsPerTFInAnnotation)
# "targets_sum_curation_effort" I actually am not sure it makes sense to use it here 
#annotationsDistributions$targets_sum_curation_effort <- annotationsDistributions$targets_sum_curation_effort / max(annotationsDistributions$targets_sum_curation_effort)

baseCols <- c("annotation","annotation_id","size")
annotationsDistributions <- fullStatsDF[,c(baseCols,annotationsDistributionsCol)]


annotationsDistributionsPivot <- pivot_longer(annotationsDistributions,cols = annotationsDistributionsCol)


ggplot(annotationsDistributionsPivot, aes(x=value, fill=name, color=name)) + geom_density(alpha=.3) + 
  facet_grid(rows = c("annotation","size"),
             scales = "free_y")
ggsave("figures_data_n_scripts/figure_annotationsDists.png", units = "cm",height = 15, width = 20, dpi = 300, bg = "white")

ggplot(annotationsDistributionsPivot, aes(x=log10(value), fill=name, color=name)) + geom_density(alpha=.3) + 
  facet_grid(rows = c("annotation","size"),
             scales = "free_y")
ggsave("figures_data_n_scripts/figure_annotationsDists_log10.png", units = "cm",height = 15, width = 20, dpi = 300, bg = "white")

## OTHERS
ggplot(annotationsDistributionsPivot, aes(x=name, y=value, fill=name, color=name)) + geom_boxplot(alpha=.3) + 
  facet_grid(rows = c("annotation","size"),
             scales = "free_y") + 
  stat_compare_means(comparisons = list(c("pTargetsW", "pTargetsF"), 
                                        c("probability", "pTargetsF"), 
                                        c("probability", "pTargetsW")),
                     size=2, method = "wilcox.test")+
  ylab("Probability | Frequency(p-value < 0.05)")+
  xlab("") +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  
ggsave("figures_data_n_scripts/figure_annotationsDists_boxplot.png", units = "cm",height = 20, width = 20, dpi = 300, bg = "white")



ggplot(annotationsDistributionsPivot,aes(x=value, y=name)) + 
  geom_violin(aes(fill=name, color=name), alpha=.3,width=1.1) +  
  geom_boxplot(aes(color=name), width=0.3, fill="grey", alpha=0.2) +
  facet_grid(rows = c("annotation","size"),
             scales = "free_y")
ggsave("figures_data_n_scripts/figure_annotationsDists_violin.png", units = "cm",height = 15, width = 20, dpi = 300, bg = "white")



