setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),".."))
getwd()
library(ggplot2)
library(gridExtra)
library(stringr)
library(smplot2)
library(ggrepel)
library(ggpubr)
library(vroom)
library(tidyverse)
library(colorspace)
library(cowplot)


swr = function(string, nwrap=50) {paste(strwrap(string, width=nwrap), collapse="\n")}; swr = Vectorize(swr)

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

dataFile <- 'random_lists_analysis/masterSimsResults.tsv'
simulationsResultsDF <- read.delim(dataFile)
simulationsResultsDF$annotation <- gsub("_"," ",simulationsResultsDF$annotation)
colnames(simulationsResultsDF)[1] <- "annotation_id"

dataFile <- 'data/observations/probabilities.csv'
annotationsSelectionProbs <- read.delim(dataFile,sep = ',')
annotationsSelectionProbs$probability

ann_info <- vroom("data/annotation_info_table.tsv")
ann_info <- rbind(ann_info,c("GO:0042493","response to xenobiotic stimulus"))

fullStatsDF <- merge(ann_info, merge(annotationsSelectionProbs,
                                     simulationsResultsDF))

simulationStats <- c("pTargetsF","pTargetsW")
probabilityStats <- c("probability")

# myX <- "rTargetsW"; myY <- "tfs_sum_curation_effort"; size <- "20"

corMethod <- "spearman"
nLabels <- 3
set.seed(9998)
myplots <- list()
mySizes <- sort(unique(fullStatsDF$size))
for (size in mySizes){
  for (myX in probabilityStats){
    for (myY in simulationStats){
      plotname <- paste0(myY,myX,size)
      print(plotname)
      subsetDF <- fullStatsDF[fullStatsDF$size == size,]
      mylabel_x <- max(subsetDF[,myX],na.rm = T) - max(subsetDF[,myX],na.rm = T) * 0.35
      mylabel_y <- max(subsetDF[,myY],na.rm = T) * 0.1
      my_ylim <- max(subsetDF[,myY],na.rm = T)
      my_xlim <- max(subsetDF[,myX],na.rm = T)

      myplots[[plotname]] <- ggplot(data = subsetDF,
                       mapping = aes(x = .data[[myX]], y = .data[[myY]],
                                     label = swr(term,30))) +
        geom_point(aes(color = annotation), shape = 20, fill='grey', stroke=0.1, size = 1) +
        sm_statCorr(color = "darkgrey", corr_method = corMethod, alternative = "two.sided",linetype = "solid", R2 = TRUE,
                    separate_by = ", ", size=1, text_size = 2.5, linewidth = 0.5, label_x = mylabel_x, label_y = mylabel_y)+
        ylim(0,min(c(1,my_ylim))+0.25) +
        xlim(0,max(c(1,my_xlim)))+
        geom_text_repel(data=subsetDF %>% arrange(desc(.data[[myY]])) %>% group_by(annotation) %>% slice_head(n=nLabels),
                        aes(x = .data[[myX]], y = .data[[myY]], color = annotation),
                        size = 2, min.segment.length=0.1, max.overlaps=50, force_pull=1, force=0.1,
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
        xlab(paste0("Size ",size)) + ylab("")
    }
  }
}

args <- c(myplots, list(ncol = 2, bottom = text_grob("Theoretical probability", size=12, face = 'bold'),
                                  top = text_grob(paste0("Fisher",paste(rep(" ",100),collapse = ""),"Wallenius"), size=12, face = 'bold'),
                                  left = text_grob("Frequency(p-value < 0.05)", size=12, face = 'bold', rot = 90)))

outName <- paste0("figures_data_n_scripts/figure_corrPermutsVSprobs.png") 
ggsave(outName, plot = do.call(grid.arrange, args), units = "cm",height = 27, width = 22, dpi = 300, bg = "white")


