setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),".."))
getwd()
library(ggplot2)
library(gridExtra)
library(stringr)
library(smplot2)
library(ggrepel)
library(ggpubr)

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
# Selected
curationEfforCol <- "targets_sum_curation_effort"

nLabels <- 20
maxProbsToLabel_prob <- annotationsCollectriStats$probability[order(annotationsCollectriStats$probability,decreasing = T)][nLabels]
selectedToLabel <- annotationsCollectriStats$probability >= maxProbsToLabel_prob

myPlot <- ggplot(data = annotationsCollectriStats, 
       mapping = aes(x = .data[[curationEfforCol]], y = probability,
                     label = annotation_id)) +
  geom_point(shape = 21, fill='grey', color = "white",stroke=0.1, size = .5) + 
  sm_statCorr(color = "black", corr_method = "pearson", linetype = "dashed", R2 = TRUE,
              separate_by = ", ", size=0.2, text_size = 2,
              label_x = mylabel_x, label_y = 0.98)+
  ylim(0,1)+
  geom_text_repel(data=subset(annotationsCollectriStats, selectedToLabel),
                  aes(x = .data[[curationEfforCol]], y = probability, color = annotation),
                  size = 1, min.segment.length=0.1, max.overlaps=20, force_pull=2, force=0.1,
                  segment.linetype = 3, segment.size=0.05) + 
  scale_color_manual(values = colors_blind)+
  theme(panel.background = element_blank(),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        legend.key.size = unit(0.2,"cm"),
        legend.key = element_rect(fill = "white"),
        plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))+
  ylab("Annotation Probability")+
  xlab("Targets Total Curation Effort")

ggsave("figure_AnotsCurEffortvsProb_selected.png", plot = myPlot, units = "cm",height = 6, width = 8,dpi = 300, bg = "white")
ggsave("figure_AnotsCurEffortvsProb_selected.tiff", plot = myPlot, units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")


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






myplot <- ggplot(annotationsCurationEffortStats, aes(x=targets_average_curation_effort, fill=annotation)) +
  geom_density() + 
  scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
  facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Average") +
  theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())

ggsave("figure_annotationsCurationEffortStats_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")

ggsave("figure_annotationsCurationEffortStats_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")


myplot <- ggplot(annotationsCurationEffortStats, aes(x=targets_median_curation_effort, fill=annotation)) +
  geom_density() + 
  scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
  facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Median") + 
  theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())

ggsave("figure_annotationsCurationEffortStats_median.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")

ggsave("figure_annotationsCurationEffortStats_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")

myplot <- ggplot(annotationsCurationEffortStats, aes(x=targets_sum_curation_effort, fill=annotation)) +
  geom_density() + 
  scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
  facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Sum") + 
  theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())

ggsave("figure_annotationsCurationEffortStats_sum.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")

ggsave("figure_annotationsCurationEffortStats_sum.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")
