library(vroom)
library(glue)
library(dplyr)
library(parallel)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(stringr)
library(data.table)
library(colorspace)
library(gtable)
library(gridExtra)
library(ggpubr)
library(sdamr)

currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))

########## Null simulations results ##########
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

annotInfo <- read.delim("data/annotation_info_table.tsv")
resultsDF <- read.delim("data/EnrResults/masterSummary.tsv")
resultsDF$annotation <- gsub("_", " ", resultsDF$annotation)
#annotInfo <- annotInfo[annotInfo$annotation_id %in% resultsDF$terms,]
fullContent <- merge(annotInfo,resultsDF,by.x ="annotation_id",by.y = "terms",all.y = T)
fullContent <- as.data.table(fullContent)

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

#### IQR Density Plots #########################################################
myplots <- list()

dark_colors <- darken(colors_blind, 0.4)
light_colors <- lighten(colors_blind,0.4)

TargetsSEAresults$annotationColors <- colors_blind[match(TargetsSEAresults$annotation,names(colors_blind))]
TargetsSEAresults$classesColors <- TargetsSEAresults$annotationColors
TargetsSEAresults$classesColors[TargetsSEAresults$typeOfAnalysis == "Target_Hypergeom"] <- darken(TargetsSEAresults$classesColors[TargetsSEAresults$typeOfAnalysis == "Target_Hypergeom"], 0.4)
TargetsSEAresults$classesColors[TargetsSEAresults$typeOfAnalysis != "Target_Hypergeom"] <- lighten(TargetsSEAresults$classesColors[TargetsSEAresults$typeOfAnalysis != "Target_Hypergeom"], 0.4)

TargetsSEAresults$classes <- paste(TargetsSEAresults$annotation,TargetsSEAresults$typeOfAnalysis)

colorClasses <- unique(TargetsSEAresults$classesColors)
names(colorClasses) <- unique(TargetsSEAresults$classes)

for (annot in annotations){
  subsetDF <- TargetsSEAresults[TargetsSEAresults$annotation == annot,] # "GO BP","KEGG","Reactome","WikiPathways"
  myplots[[annot]] <- ggplot(subsetDF, aes(x=RankingIQR, color=classes)) + 
                      geom_density(color=subsetDF$classesColors, fill=subsetDF$classesColors) + 
                      xlab(annot) + ylab("")
}
args <- c(myplots, list(ncol = 2, bottom="Rank IQR",left="Density"))
ggsave("figure_IQRdist.tiff", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 500, compression = "lzw",bg = "white")
ggsave("figure_IQRdist.png", plot = do.call(grid.arrange, args),units = "cm",height = 12, width = 18,dpi = 400, bg = "white")

################################################################################
################################################################################

# TO EXPLORE ### IQR in Range Points #################################################

annot = "KEGG"
subsetDF <- TargetsSEAresults[TargetsSEAresults$annotation == annot,]
subsetDF$RankingIQRscaled <- subsetDF$RankingIQR / nrow(subsetDF) * 100
subsetDF$ymin <- subsetDF$RankingMedian - subsetDF$RankingIQRscaled*0.5
subsetDF$ymin[subsetDF$ymin <= 1] <- 1
subsetDF$ymax <- subsetDF$RankingMedian + subsetDF$RankingIQRscaled*0.5
ggplot(subsetDF, aes(x = size, y = RankingMedian)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax), 
                  position=position_jitter(), 
                  linetype='dotted') 

################################################################################
################################################################################


#### Histogram of Results ######################################################

fullContent <- TargetsSEAresults
fullContent <- fullContent[!is.na(fullContent$p_times),]
mysizes <- sort(unique(fullContent$size))
annotation = "GO BP"

fullPlots <- lapply(annotations, function(annotation){
  base_color <- colors_blind[names(colors_blind) %in% annotation]
  dark_color <- darken(base_color, 0.4)
  light_color <- lighten(base_color,0.4)
  color_vals_insert <- c(dark_color, light_color)
  names(color_vals_insert) <- c("Fisher's Exact Test","Wallenius Test")
  
  fullContentAnn <- fullContent %>% filter(annotation == !!annotation, size %in% mysizes) %>%
    arrange(size) %>%
    mutate(typeOfAnalysis = ifelse(typeOfAnalysis == "Target_Hypergeom","Fisher's Exact Test","Wallenius Test"),
           typeOfAnalysis = factor(typeOfAnalysis, levels = c("Fisher's Exact Test","Wallenius Test")),
           sizeName = factor(paste(size, "TFs"), levels = unique(paste(size, "TFs"))),
           percentage = p_times * 100)
  
  hist_info <- rbindlist(lapply(levels(fullContentAnn$sizeName), function(size){
    hist_info <- fullContentAnn %>% filter(sizeName == !!size)
    hist_info <- rbindlist(lapply(levels(hist_info$typeOfAnalysis), function(typeOfAnalysis){
      hist_info <- hist_info %>% filter(typeOfAnalysis == !!typeOfAnalysis)
      p <- ggplot(hist_info, aes(x = percentage))+
        geom_histogram(breaks = seq(0,100,by = 2))
      hist_info <- ggplot_build(p)$data[[1]]
      data.frame(counts = hist_info$count,percentage = hist_info$x, xmin = hist_info$xmin, xmax = hist_info$xmax,
                 typeOfAnalysis = typeOfAnalysis, sizeName = size) %>% filter(counts != 0)
    }))
  })) %>% mutate(typeOfAnalysis = factor(typeOfAnalysis, levels = c("Fisher's Exact Test","Wallenius Test")),
                 sizeName = factor(sizeName, levels = levels(fullContentAnn$sizeName)))
  
  p <- ggplot(hist_info)+
    facet_wrap(~sizeName, ncol = 2, scales = "free")+
    geom_col(aes(x = percentage, y = counts, fill = typeOfAnalysis),linewidth = 0.2, position = "identity", color = "black", alpha = 0.6, width = 2,stroke=0)+
    #geom_vline(xintercept = 0.05, linewidth = 0.2, color = "red", linetype = "dotted")+
    scale_fill_manual(values = color_vals_insert)+
    theme(panel.background = element_blank(),
          strip.text = element_text(size = 4),
          strip.background = element_rect(fill = "white"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 4),
          title = element_text(face = "bold"),
          axis.title = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 4),
          legend.key.size = unit(0.2,"cm"))+
    ylab("Frequency")+
    xlab("Percentage of times significant")+
    scale_y_log10()+
    labs(fill = annotation)+
    scale_x_continuous(limits = c(0,100))
    
  #p <- shift_legend2(p)
  #dev.off()
  
  FisherResults <- fullContentAnn[fullContentAnn$typeOfAnalysis == "Fisher's Exact Test",]
  FisherResults <- FisherResults[order(FisherResults$RankingMedian),]
  topAnnotationsFisher <- unique(FisherResults$annotation_id)[1:10]
  content_top <- fullContentAnn %>% filter(annotation_id %in% topAnnotationsFisher)
  
  content_top$ymin <- content_top$RankingMedian - content_top$RankingMedian*0.5
  content_top$ymin[content_top$ymin <= 1] <- 1
  content_top$ymax <- content_top$RankingMedian + content_top$RankingMedian*0.5
  
  content_top$size <- factor(content_top$size,levels = sort(unique(content_top$size)))

  p2 <- ggplot(content_top, aes(x = size, y = RankingMedian))+
    #geom_jitter(aes(color = typeOfAnalysis), size = 0.3, width = 0.2)+
    geom_pointrange(aes(color = typeOfAnalysis, ymin = ymin, ymax = ymax), 
                    position=position_jitter(), linetype='solid',
                    size = 0.001, fatten = 0.01) +
    scale_color_manual(values = color_vals_insert)+
    theme(panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 4),
          title = element_text(face = "bold"),
          strip.text = element_text(size = 4,face = "bold"),
          axis.title = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          legend.margin=margin(0,8,0,0),
          legend.box.margin=margin(-10,0,-10,-10),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          legend.position = "none",
          plot.title = element_text(size = 4, hjust = 0.5))+
    ggtitle(annotation)+
    xlab("Number of TFs")+
    ylab("Median Rank of Top10 Terms")
  
  return(list(p,p2))
})

gg1 <- plot_grid(plotlist = sapply(fullPlots,"[",1))
gg2 <- plot_grid(plotlist = sapply(fullPlots,"[",2), nrow = 1)

gg4_v2 <- ggdraw()+
  draw_plot(gg1,x = 0, y = 0.2, height = 0.8, width = 1)+
  draw_plot(gg2, x = 0, y = 0, height = 0.2, width = 1) +
  draw_plot_label(c("A","B"), x = c(0, 0), y = c(1, 0.225), fontface = "plain", family = "serif", size = 10)

# gg4_v2 <- plot_grid(gg4_v2, gg_legend, rel_widths = c(0.9,0.2))
ggsave("MegaFigureTargets.tiff", plot = gg4_v2, units = "cm",height = 18, width = 20,dpi = 600, compression = "lzw", bg = "white")
ggsave("MegaFigureTargets.png", plot = gg4_v2, units = "cm",height = 18, width = 20,dpi = 400, bg = "white")


