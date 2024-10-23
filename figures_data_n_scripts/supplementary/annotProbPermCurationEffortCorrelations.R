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

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

dataFile <- 'data/observations/probabilities.csv'
annotationsSelectionProbs <- read.delim(dataFile,sep = ',')
annotationsSelectionProbs$probability

dataFile <- 'random_lists_analysis/masterSimsResults.tsv'
simulationsResultsDF <- read.delim(dataFile)
simulationsResultsDF$annotation <- gsub("_"," ",simulationsResultsDF$annotation)
colnames(simulationsResultsDF)[1] <- "annotation_id"

dataFile <- 'data/observations/annotationsCurationEffortStats.tsv'
annotationsCurationEffortStats <- read.delim(dataFile)
annotationsCurationEffortStats$annotation <- gsub("_"," ",annotationsCurationEffortStats$annotation)
nrow(annotationsCurationEffortStats)

ann_info <- vroom("data/annotation_info_table.tsv")
ann_info <- rbind(ann_info,c("GO:0042493","response to xenobiotic stimulus"))

# FIRST CHECK --> ok
# table(annotationsCurationEffortStats$annotation) == table(annotationsSelectionProbs$ann)
# all(annotationsCurationEffortStats$annotation_id %in% annotationsSelectionProbs$annotation_id)
# all(annotationsSelectionProbs$annotation_id %in% annotationsCurationEffortStats$annotation_id)

# SECOND CHECK --> ERROR
# table(annotationsCurationEffortStats$annotation) == table(simulationsResultsDF$annotation) / length(unique(simulationsResultsDF$size))


##  THERE ARE MORE ANNOTATIONS IN THE SIMULATIONS RESULTS THAN IN THE OTHER RESULTS 
##  WHAT ARE WE MISISNG or DUPLICATING?
## any(table(simulationsResultsDF$annotation_id) > length(unique(simulationsResultsDF$size)))
## NOPE! IS NOT DUPLICATING
# TRUE  # all(annotationsCurationEffortStats$annotation_id %in% simulationsResultsDF$annotation_id)
# FALSE # all(simulationsResultsDF$annotation_id %in% annotationsCurationEffortStats$annotation_id)
## SIMULATION RESULTS INCLUDES TERM
## unique(setdiff(simulationsResultsDF$annotation_id,annotationsCurationEffortStats$annotation_id))
missingIDs <- unique(setdiff(simulationsResultsDF$annotation_id,annotationsCurationEffortStats$annotation_id))
View(simulationsResultsDF[simulationsResultsDF$annotation_id %in% missingIDs,])

# TRUE # all(is.na(simulationsResultsDF[simulationsResultsDF$annotation_id %in% missingIDs,]$pTFs))
# FALSE # all(is.na(simulationsResultsDF[simulationsResultsDF$annotation_id %in% missingIDs,]$pTargetsF))
## FOR SOME REASON THERE ARE RESULTS IN TFs SEA NOT ASSOCIATED WITH PROBABILITY/CURATIONEFFORT
## IS IT BECAUSE WE ONLY USED THE UNIVERSE OF TARGETS TO CALCULATE PROBABILITY/CURATIONEFFORT ???

fullStatsDF <- merge(ann_info, 
                     merge(annotationsCurationEffortStats,
                           merge(annotationsSelectionProbs,simulationsResultsDF)
                           )
                     )

totalAnnotationPerDB <- table(simulationsResultsDF[simulationsResultsDF$size == "10",]$annotation)

colnames(annotationsSelectionProbs); colnames(simulationsResultsDF); colnames(annotationsCurationEffortStats)

table(simulationsResultsDF$size)
table(fullStatsDF$size)

simulationStats <- c("pTFs","rTFs","pTargetsF","rTargetsF","pTargetsW","rTargetsW")
curationStats <- c("targets_sum_curation_effort","tfs_sum_curation_effort")

correlationVariables <- list(list(Y="probability",X=simulationStats),
                             list(Y="probability",X=curationStats),
                             list(Y="targets_sum_curation_effort",X=simulationStats),
                             list(Y="tfs_sum_curation_effort",X=simulationStats))

# myX <- "rTargetsW"; myY <- "tfs_sum_curation_effort"; size <- "20"

nLabels <- 20
set.seed(9998)
rankingRelative <- T
correlationStatistics <- c("file","myY","myX","size","correlation_coefficient","R2","p_value")
for (comparison in 1:length(correlationVariables)){
  myY <- correlationVariables[[comparison]][["Y"]]
  for (myX in correlationVariables[[comparison]][["X"]]){
    mySizes <- unique(fullStatsDF$size)
    if (grepl("effort",myX)){
      mySizes <- c("3") # ÑAPA TO IGNORE SIZE WHEN CORRELATING probability ~ curationStats, SIZE ONLY AFFECTS simulationStats
    }else{
      mySizes <- unique(fullStatsDF$size)
    }
    for (size in mySizes){
      print(paste(myY,myX,size))
      
      subsetDF <- fullStatsDF[fullStatsDF$size == size,]
      
      mylabel_x <- max(subsetDF[,myX],na.rm = T) - max(subsetDF[,myX],na.rm = T) * 0.2
      mylabel_y <- max(subsetDF[,myY],na.rm = T) * 0.5
      my_ylim <- max(subsetDF[,myY],na.rm = T)
      my_xlim <- max(subsetDF[,myX],na.rm = T)
      tag <- ""
      if (myY == "probability"){
        maxProbsToLabel <- sort(subsetDF[,myY],decreasing = T)[nLabels]
        selectedToLabel <- subsetDF[,myY] >= maxProbsToLabel
      }else if (grepl("rT",myX)){
        maxProbsToLabel <- sort(subsetDF[,myX],decreasing = F)[nLabels]
        selectedToLabel <- subsetDF[,myX] <= maxProbsToLabel
        if(rankingRelative){
          my_xlim <- 100
          tag <- "relative"
          for (myannotation in unique(subsetDF$annotation)){
            subsetDF[,myX][subsetDF$annotation == myannotation] <- subsetDF[,myX][subsetDF$annotation == myannotation] / totalAnnotationPerDB[[myannotation]] * 100
          }
        }
      }else{
        maxProbsToLabel <- sort(subsetDF[,myY],decreasing = T)[nLabels]
        selectedToLabel <- subsetDF[,myY] >= maxProbsToLabel
      }
        
      # Calculate correlation and extract statistics
      corMethod <- "spearman"
      cor_result <- cor.test(subsetDF[,myX], subsetDF[,myY], method = corMethod,
                             alternative = "two.sided")
      
      # Extracting statistics
      correlation_coefficient <- cor_result$estimate # Correlation coefficient
      R2 <- cor_result$estimate ^ 2
      p_value <- cor_result$p.value                      # p-value
      
      myPlot <- ggplot(data = subsetDF,
                       mapping = aes(x = .data[[myX]], y = .data[[myY]],
                                     label = swr(term,30))) +
        geom_point(aes(color = annotation), shape = 20, fill='grey', stroke=0.1, size = 1) +
        sm_statCorr(color = "black", corr_method = corMethod, alternative = "two.sided",linetype = "dashed", R2 = TRUE,
                    separate_by = ", ", size=1, text_size = 5,
                    label_x = mylabel_x, label_y = mylabel_y)+
        ylim(0,max(c(1,my_ylim)))+
        xlim(0,max(c(1,my_xlim)))+
        geom_text_repel(data=subset(subsetDF, selectedToLabel),
                        aes(x = .data[[myX]], y = .data[[myY]], color = annotation),
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
        xlab(paste0(myX," at size ",size))+
        ylab(myY)
      
      outName <- paste0("figures_data_n_scripts/annotProbPermCurationEffortCorrelations/figure_correlation",tag,myX,myY,"_size",size,".png") 
      correlationStatistics <- rbind(correlationStatistics,cbind(basename(outName),myY,myX,size,correlation_coefficient,R2,p_value))
      ggsave(outName, plot = myPlot, units = "cm",height = 20, width = 30, dpi = 300, bg = "white")
    }
  }
}
correlationStatistics <- as.data.frame(correlationStatistics)
# colnames(correlationStatistics) <- c("file","myY","myX","size","correlation_coefficient","p_value","confidence_interval")

write.table(correlationStatistics,"figures_data_n_scripts/annotProbPermCurationEffortCorrelations/correlation_table.tsv",sep = "\t",quote = F,row.names = F,col.names = F)



# Selected

simulationStats <- c("pTargetsF","pTargetsW")
curationStats <- c("targets_sum_curation_effort")

correlationVariables <- list(list(Y="probability",X=simulationStats),
                             list(Y="probability",X=curationStats),
                             list(Y="targets_sum_curation_effort",X=simulationStats),
                             list(Y="tfs_sum_curation_effort",X=simulationStats))

# myX <- "rTargetsW"; myY <- "tfs_sum_curation_effort"; size <- "20"

nLabels <- 20
set.seed(9998)
rankingRelative <- T
correlationStatistics <- c("file","myY","myX","size","correlation_coefficient","R2","p_value")
for (comparison in 1:length(correlationVariables)){
  myY <- correlationVariables[[comparison]][["Y"]]
  for (myX in correlationVariables[[comparison]][["X"]]){
    mySizes <- unique(fullStatsDF$size)
    if (grepl("effort",myX)){
      mySizes <- c("3") # ÑAPA TO IGNORE SIZE WHEN CORRELATING probability ~ curationStats, SIZE ONLY AFFECTS simulationStats
    }else{
      mySizes <- unique(fullStatsDF$size)
    }
    for (size in mySizes){
      print(paste(myY,myX,size))
      
      subsetDF <- fullStatsDF[fullStatsDF$size == size,]
      
      mylabel_x <- max(subsetDF[,myX],na.rm = T) - max(subsetDF[,myX],na.rm = T) * 0.2
      mylabel_y <- max(subsetDF[,myY],na.rm = T) * 0.5
      my_ylim <- max(subsetDF[,myY],na.rm = T)
      my_xlim <- max(subsetDF[,myX],na.rm = T)
      tag <- ""
      if (myY == "probability"){
        maxProbsToLabel <- sort(subsetDF[,myY],decreasing = T)[nLabels]
        selectedToLabel <- subsetDF[,myY] >= maxProbsToLabel
      }else if (grepl("rT",myX)){
        maxProbsToLabel <- sort(subsetDF[,myX],decreasing = F)[nLabels]
        selectedToLabel <- subsetDF[,myX] <= maxProbsToLabel
        if(rankingRelative){
          my_xlim <- 100
          tag <- "relative"
          for (myannotation in unique(subsetDF$annotation)){
            subsetDF[,myX][subsetDF$annotation == myannotation] <- subsetDF[,myX][subsetDF$annotation == myannotation] / totalAnnotationPerDB[[myannotation]] * 100
          }
        }
      }else{
        maxProbsToLabel <- sort(subsetDF[,myY],decreasing = T)[nLabels]
        selectedToLabel <- subsetDF[,myY] >= maxProbsToLabel
      }
      
      # Calculate correlation and extract statistics
      corMethod <- "spearman"
      cor_result <- cor.test(subsetDF[,myX], subsetDF[,myY], method = corMethod,
                             alternative = "two.sided")
      
      # Extracting statistics
      correlation_coefficient <- cor_result$estimate # Correlation coefficient
      R2 <- cor_result$estimate ^ 2
      p_value <- cor_result$p.value                      # p-value
      
      myPlot <- ggplot(data = subsetDF,
                       mapping = aes(x = .data[[myX]], y = .data[[myY]],
                                     label = swr(term,30))) +
        geom_point(aes(color = annotation), shape = 20, fill='grey', stroke=0.1, size = 1) +
        sm_statCorr(color = "black", corr_method = corMethod, alternative = "two.sided",linetype = "dashed", R2 = TRUE,
                    separate_by = ", ", size=1, text_size = 5,
                    label_x = mylabel_x, label_y = mylabel_y)+
        ylim(0,max(c(1,my_ylim)))+
        xlim(0,max(c(1,my_xlim)))+
        geom_text_repel(data=subset(subsetDF, selectedToLabel),
                        aes(x = .data[[myX]], y = .data[[myY]], color = annotation),
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
        xlab(paste0(myX," at size ",size))+
        ylab(myY)
      
      outName <- paste0("figures_data_n_scripts/annotProbPermCurationEffortCorrelations/figure_correlation",tag,myX,myY,"_size",size,".png") 
      correlationStatistics <- rbind(correlationStatistics,cbind(basename(outName),myY,myX,size,correlation_coefficient,R2,p_value))
      ggsave(outName, plot = myPlot, units = "cm",height = 20, width = 30, dpi = 300, bg = "white")
    }
  }
}





# toCorrelateStats <- colnames(annotationsCollectriStats)[grep("effort",colnames(annotationsCollectriStats))]
#minProbsToLabel_prob <- annotationsCollectriStats$probability[order(annotationsCollectriStats$probability,decreasing = F)][nLabels]
# myPlots <- list()
# for (toCorrelateStat in toCorrelateStats){
  #maxProbsToLabel_cureff <- annotationsCollectriStats[,toCorrelateStat][order(annotationsCollectriStats[,toCorrelateStat],decreasing = T)][nLabels]
  # selectedToLabel <- annotationsCollectriStats$probability >= maxProbsToLabel_prob # | 
                     # annotationsCollectriStats$probability <= minProbsToLabel_prob | 
                     # annotationsCollectriStats[,toCorrelateStat] >= maxProbsToLabel_cureff
#   
#   length(annotationsCollectriStats$annotation_id[selectedToLabel])
#   
#   mylabel_x <- min(annotationsCollectriStats[,toCorrelateStat],na.rm = T)
#   myPlots[[toCorrelateStat]] <- ggplot(data = annotationsCollectriStats, 
#                                         mapping = aes(x = .data[[toCorrelateStat]], y = probability,
#                                                       label = annotation_id)) +
#     geom_point(shape = 21, fill='grey', color = "white",stroke=0.1, size = .5) + 
#     sm_statCorr(color = "black", corr_method = "pearson", linetype = "solid", R2 = TRUE,
#                 separate_by = ", ", size=0.2, text_size = 1,
#                 label_x = mylabel_x, label_y = 0.99)+
#     ylim(0,1)+
#     geom_text_repel(data=subset(annotationsCollectriStats, selectedToLabel),
#               aes(x = .data[[toCorrelateStat]], y = probability, color = annotation),
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
#     xlab(str_to_title(gsub("_"," ",toCorrelateStat))) 
# }
#topAnnot <- text_grob("Anotations Curation Effort vs Probability", size=6, face = 'bold')
#lefAnnot <- text_grob("Annotation Probability", size=6, face = 'bold',rot = 90)
# args <- c(myPlots, list(ncol = 3, top="Anotations Curation Effort vs Probability",left="Annotation Probability"))
# ggsave("figure_AnotsCurEffortvsProb_latest.tiff", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 500, compression = "lzw",bg = "white")
# ggsave("figure_AnotsCurEffortvsProb_latest.png", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 400, bg = "white")


###
### AN OLD VERSION OF PREVIOUS PLOT BUT SEPARATING BY ANNOTATION DATABASE - DISCARDED AT THE MOMENT
# myPlots <- list()
# for (toCorrelateStat in toCorrelateStats){
#   myPlots[[toCorrelateStat]] <- ggplot(data = annotationsCollectriStats, 
#                                         mapping = aes(x = .data[[toCorrelateStat]], y = probability,
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
#     xlab(str_to_title(gsub("_"," ",toCorrelateStat))) 
# }
# 
# args <- c(myPlots, list(ncol = 3, top="Anotations Curation Effort vs Probability",left="Annotation Probability"))
# ggsave("figure_AnotsCurEffortvsProb_BYdatabase.tiff", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 500, compression = "lzw",bg = "white")
# ggsave("figure_AnotsCurEffortvsProb_BYdatabase.png", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 400, bg = "white")

 
# myplot <- ggplot(annotationsSelectionProbs, aes(x=targets_average_curation_effort, fill=annotation)) +
#   geom_density() + 
#   scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
#   facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Average") +
#   theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())
# 
# ggsave("figure_annotationsSelectionProbs_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")
# 
# ggsave("figure_annotationsSelectionProbs_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")
# 
# 
# myplot <- ggplot(annotationsSelectionProbs, aes(x=targets_median_curation_effort, fill=annotation)) +
#   geom_density() + 
#   scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
#   facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Median") + 
#   theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())
# 
# ggsave("figure_annotationsSelectionProbs_median.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")
# 
# ggsave("figure_annotationsSelectionProbs_average.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")
# 
# myplot <- ggplot(annotationsSelectionProbs, aes(x=targets_sum_curation_effort, fill=annotation)) +
#   geom_density() + 
#   scale_fill_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') + 
#   facet_grid(rows = vars(annotation)) + ylab("Annotation Density") + xlab("Targets Curation Effort Sum") + 
#   theme(panel.background = element_blank(), strip.background = element_blank(), strip.text = element_blank())
# 
# ggsave("figure_annotationsSelectionProbs_sum.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")
# 
# ggsave("figure_annotationsSelectionProbs_sum.tiff", plot = myplot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")
