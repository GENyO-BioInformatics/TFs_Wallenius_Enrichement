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
swr = function(string, nwrap=50) {paste(strwrap(string, width=nwrap), collapse="\n")}; swr = Vectorize(swr)


################################################################################
### 1 GENERATE DATA FOR PLOT
################################################################################

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
# collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
# annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv'); annotationFile <- annotationFiles[1]
# topKnown <- 50
# for (annotationFile in annotationFiles){
#   annotation <- gsub('.tsv','',basename(annotationFile))
#   annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
#   annotationDF <- annotationDF[annotationDF$organism == 9606,]
#   cat("Annotation:", annotation, "\n")
#   for (top in (topKnown-20):topKnown){
#     curationEffortCutOff <- unique(collectriTFsGRN$curation_effort[order(collectriTFsGRN$curation_effort,decreasing = T)])[top]
#     mostKnownTargets <- unique(collectriTFsGRN$target_genesymbol[collectriTFsGRN$curation_effort >= curationEffortCutOff])
#     mostKnownTFs <- unique(collectriTFsGRN$source_genesymbol[collectriTFsGRN$curation_effort >= curationEffortCutOff])
#     annotationsCurationEffortStats <- c()
#     # cat("Top:", top,
#     #     "\tTarget annotations", 
#     #     length(unique(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTargets])),
#     #     "\tand\tTFs annotations", 
#     #     length(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTFs]),
#     #     "\n"
#     # )
#     cat("Top: ", top,
#         "\tTarget most common annotation: ",
#         max(table(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTargets])),
#         "/",
#         mean(table(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTargets])),
#         "\tTFs most common annotation: ",
#         max(table(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTFs])),
#         "/",
#         mean(table(annotationDF$annotation_id[annotationDF$symbol %in% mostKnownTFs])),
#         "\n")
#   }
# }

################################################################################
### 2. PLOT DATA
################################################################################

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

dataFile <- 'data/observations/annotationsCurationEffortStats.tsv'
annotationsCurationEffortStats <- read.delim(dataFile)
annotationsCurationEffortStats$annotation <- gsub("_"," ",annotationsCurationEffortStats$annotation)
nrow(annotationsCurationEffortStats)

dataFile <- 'data/observations/probabilities.csv'
annotationsSelectionProbs <- read.delim(dataFile,sep = ',')

if(nrow(annotationsCurationEffortStats) != nrow(annotationsSelectionProbs)){
  cat("WTF!!!! NO PUEDE SER. ¿COMO ESTAMOS SACANDO DIFERENTE N DE ANNOTS?")
}

annotationsCollectriStats <- merge(annotationsCurationEffortStats,annotationsSelectionProbs)
ann_info <- vroom("data/annotation_info_table.tsv")

annotationsCollectriStats <- merge(ann_info,annotationsCollectriStats,all.y = T)

#ann_info <- rbind(ann_info,c("GO:0042493","response to xenobiotic stimulus"))


View(annotationsCollectriStats)

# Selected
curationEfforCol <- "targets_sum_curation_effort"

mylabel_x <- min(annotationsCollectriStats[,curationEfforCol],na.rm = T)

nLabels <- 20
maxProbsToLabel_prob <- annotationsCollectriStats$probability[order(annotationsCollectriStats$probability,decreasing = T)][nLabels]
selectedToLabel <- annotationsCollectriStats$probability >= maxProbsToLabel_prob
set.seed(999)
myPlot <- ggplot(data = annotationsCollectriStats, 
       mapping = aes(x = .data[[curationEfforCol]], y = probability,
                     label = swr(term,30))) +
  geom_point(aes(color = annotation), shape = 20, fill='grey', stroke=0.1, size = 1) + 
  sm_statCorr(color = "black", corr_method = "pearson", linetype = "dashed", R2 = TRUE,
              separate_by = ", ", size=0.2, text_size = 6,
              label_x = mylabel_x, label_y = 0.98)+
  ylim(0,1)+
  geom_text_repel(data=subset(annotationsCollectriStats, selectedToLabel),
                  aes(x = .data[[curationEfforCol]], y = probability, color = annotation),
                  size = 3.5, min.segment.length=0.1, max.overlaps=20, force_pull=2, force=0.1,
                  segment.linetype = 3, segment.size=0.05) + 
  scale_color_manual(values = colors_blind)+
  theme(panel.background = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        legend.key.size = unit(0.2,"cm"),
        legend.key = element_rect(fill = "white"),
        plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))+
  ylab("Annotation Probability")+
  xlab("Targets Total Curation Effort")

ggsave("figures_data_n_scripts/figure_AnotsCurEffortvsProb.png", plot = myPlot, units = "cm",height = 20, width = 30, dpi = 300, bg = "white")
ggsave("figures_data_n_scripts/figure_AnotsCurEffortvsProb.tiff", plot = myPlot, units = "cm",height = 15, width = 15,dpi = 500, compression = "lzw",bg = "white")


# curationEfforCols <- colnames(annotationsCollectriStats)[grep("effort",colnames(annotationsCollectriStats))]
#minProbsToLabel_prob <- annotationsCollectriStats$probability[order(annotationsCollectriStats$probability,decreasing = F)][nLabels]
# myPlots <- list()
# for (curationEfforCol in curationEfforCols){
  #maxProbsToLabel_cureff <- annotationsCollectriStats[,curationEfforCol][order(annotationsCollectriStats[,curationEfforCol],decreasing = T)][nLabels]
  # selectedToLabel <- annotationsCollectriStats$probability >= maxProbsToLabel_prob # | 
                     # annotationsCollectriStats$probability <= minProbsToLabel_prob | 
                     # annotationsCollectriStats[,curationEfforCol] >= maxProbsToLabel_cureff
#   
#   length(annotationsCollectriStats$annotation_id[selectedToLabel])
#   
#   mylabel_x <- min(annotationsCollectriStats[,curationEfforCol],na.rm = T)
#   myPlots[[curationEfforCol]] <- ggplot(data = annotationsCollectriStats, 
#                                         mapping = aes(x = .data[[curationEfforCol]], y = probability,
#                                                       label = annotation_id)) +
#     geom_point(shape = 21, fill='grey', color = "white",stroke=0.1, size = .5) + 
#     sm_statCorr(color = "black", corr_method = "pearson", linetype = "solid", R2 = TRUE,
#                 separate_by = ", ", size=0.2, text_size = 1,
#                 label_x = mylabel_x, label_y = 0.99)+
#     ylim(0,1)+
#     geom_text_repel(data=subset(annotationsCollectriStats, selectedToLabel),
#               aes(x = .data[[curationEfforCol]], y = probability, color = annotation),
#               size = 0.8, min.segment.length=0.1, max.overlaps=20, force_pull=2, force=0.1,
#               segment.linetype = 3, segment.size=0.05) + 
#     scale_color_manual(values = colors_blind)+
#     theme(panel.background = element_blank(),
#           axis.title = element_text(size = 5),
#           axis.text = element_text(size = 4),
#           legend.title = element_text(size = 4),
#           legend.text = element_text(size = 4),
#           legend.key.size = unit(0.2,"cm"),
#           legend.key = element_rect(fill = "white"),
#           plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))+
#     ylab("")+
#     xlab(str_to_title(gsub("_"," ",curationEfforCol))) 
# }
#topAnnot <- text_grob("Anotations Curation Effort vs Probability", size=6, face = 'bold')
#lefAnnot <- text_grob("Annotation Probability", size=6, face = 'bold',rot = 90)
# args <- c(myPlots, list(ncol = 3, top="Anotations Curation Effort vs Probability",left="Annotation Probability"))
# ggsave("figure_AnotsCurEffortvsProb_latest.tiff", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 500, compression = "lzw",bg = "white")
# ggsave("figure_AnotsCurEffortvsProb_latest.png", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 400, bg = "white")


###
### AN OLD VERSION OF PREVIOUS PLOT BUT SEPARATING BY ANNOTATION DATABASE - DISCARDED AT THE MOMENT
# myPlots <- list()
# for (curationEfforCol in curationEfforCols){
#   myPlots[[curationEfforCol]] <- ggplot(data = annotationsCollectriStats, 
#                                         mapping = aes(x = .data[[curationEfforCol]], y = probability,
#                                                       label = annotation_id, fill = annotation)) +
#     geom_point(shape = 21, color = "white",stroke=0.1, size = .5) + 
#     scale_fill_manual(values = colors_blind) +
#     sm_statCorr(color = "black", corr_method = "pearson", linetype = "solid", R2 = TRUE,
#                 separate_by = "\n", size=0.2,text_size = 2)+
#     theme(panel.background = element_blank(),
#           axis.title = element_text(size = 5),
#           axis.text = element_text(size = 4),
#           legend.title = element_text(size = 4),
#           legend.text = element_text(size = 4),
#           legend.key.size = unit(0.2,"cm"),
#           legend.key = element_rect(fill = "white"),
#           plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))+
#     ylab("")+
#     xlab(str_to_title(gsub("_"," ",curationEfforCol))) 
# }
# 
# args <- c(myPlots, list(ncol = 3, top="Anotations Curation Effort vs Probability",left="Annotation Probability"))
# ggsave("figure_AnotsCurEffortvsProb_BYdatabase.tiff", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 500, compression = "lzw",bg = "white")
# ggsave("figure_AnotsCurEffortvsProb_BYdatabase.png", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 400, bg = "white")

 
# myplot <- ggplot(annotationsCurationEffortStats, aes(x=targets_average_curation_effort, fill=annotation)) +
#   geom_density() + 
#   scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
#   facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Average") +
#   theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())
# 
# ggsave("figure_annotationsCurationEffortStats_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")
# 
# ggsave("figure_annotationsCurationEffortStats_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")
# 
# 
# myplot <- ggplot(annotationsCurationEffortStats, aes(x=targets_median_curation_effort, fill=annotation)) +
#   geom_density() + 
#   scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
#   facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Median") + 
#   theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())
# 
# ggsave("figure_annotationsCurationEffortStats_median.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")
# 
# ggsave("figure_annotationsCurationEffortStats_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")
# 
# myplot <- ggplot(annotationsCurationEffortStats, aes(x=targets_sum_curation_effort, fill=annotation)) +
#   geom_density() + 
#   scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
#   facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Sum") + 
#   theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())
# 
# ggsave("figure_annotationsCurationEffortStats_sum.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")
# 
# ggsave("figure_annotationsCurationEffortStats_sum.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")
