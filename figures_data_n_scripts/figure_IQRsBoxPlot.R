library(vroom)
library(ggplot2)
library(ggrepel)
#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(tidyverse)
library(colorspace)
library(ggrepel)
library(cowplot)
set.seed(1234)

currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")
l_color <- lighten(colors_blind, 0.2)
d_color <- darken(colors_blind, 0.2)
color_pal <- c("Hypergeometric" = d_color, "Non-Central Hypergeometric" = l_color)


TargetsSEA <- read.delim("data/EnrResults/TargetsSEAresultsMelted.tsv")
TargetsSEA <- TargetsSEA[grep("Target_",TargetsSEA$typeOfAnalysis),] 
TargetsSEA$annotation <- gsub('_',' ',TargetsSEA$annotation)
TargetsSEA$RankingMedian
table(TargetsSEA$typeOfAnalysis)


topToShow <- 15

topBiassedAnnotsDF <- TargetsSEA[TargetsSEA$typeOfAnalysis == "Target_Hypergeom",] %>% arrange(RankingMedian) %>% group_by(annotation,size) %>% slice_head(n=topToShow)
topBiassedAnnotsDF <- TargetsSEA[TargetsSEA$annotation_id %in% unique(topBiassedAnnotsDF$annotation_id),]

topBiassedAnnotsDF$IQR <- topBiassedAnnotsDF$RankingQ3 - topBiassedAnnotsDF$RankingQ1

write.table(topBiassedAnnotsDF,file = "figures_data_n_scripts/RankingPlots_topBiassedAnnots.tsv",sep = "\t",quote = F,row.names = F,col.names = T)

topBiassedAnnotsDF <- topBiassedAnnotsDF %>% mutate(typeOfAnalysis = ifelse(typeOfAnalysis == "Target_Hypergeom", "Hypergeometric", "Non-Central Hypergeometric"))

topBiassedAnnotsDF$class <- paste0(topBiassedAnnotsDF$typeOfAnalysis,".",topBiassedAnnotsDF$annotation)
topBiassedAnnotsDF$colors <- color_pal[match(topBiassedAnnotsDF$class,names(color_pal))]

topBiassedAnnotsDF <- topBiassedAnnotsDF %>% mutate(annotation = factor(annotation, levels = names(colors_blind)), size=factor(size,levels = sort(unique(size)))) 

topBiassedAnnotsDF <- topBiassedAnnotsDF %>% group_by(typeOfAnalysis) %>% mutate(n = n(), x_vals = runif(n, -0.2, 0.2)) %>% ungroup() %>%  mutate(x_vals = ifelse(typeOfAnalysis == "Hypergeometric", x_vals -0.25, x_vals + 0.25))


annotRankingTops <- list()
for (ann in levels(topBiassedAnnotsDF$annotation)){
  sub_topBiassedAnnotsDF <- topBiassedAnnotsDF[topBiassedAnnotsDF$annotation == ann,]
  annotRankingTops[[ann]] <- ggplot(data = sub_topBiassedAnnotsDF, aes(y = RankingMedian)) +
      geom_point(aes(x = as.numeric(size) + x_vals, y = RankingMedian,
                     color = class, fill = class),
                     size = 0.5, stroke = 0) + 
      geom_pointrange(aes(x = as.numeric(factor(size)) + x_vals,
                        color = class, fill = class, ymin = RankingQ1, ymax = RankingQ3),
                    linetype='solid',
                    size = 0.001, fatten = 0.01) +
    scale_x_discrete(labels = levels(sub_topBiassedAnnotsDF$size))+
    theme_classic() +
    theme(legend.position = "bottom", text = element_text(size = 6)) +
    scale_fill_manual(values = color_pal,guide = 'none') +
    scale_color_manual(values = color_pal,guide = 'none') + 
    geom_text_repel(data = sub_topBiassedAnnotsDF,
    aes(x = as.numeric(factor(size)) + x_vals,
        y = RankingMedian,
        label = swr(term,15)),
    seed = 343435,
    min.segment.length = unit(0, "lines"),
    segment.size = 0.2,
    size = 2,
    lineheight = 0.6,
    box.padding = 0.5,
    force_pull = 1) + 
    xlab("Size of TFs Lists") +
    ylab("Median Ranking Position") 
}
gg <- plot_grid(plotlist = annotRankingTops)
ggsave("figures_data_n_scripts/RankingRandomTFtargetsLists_terms.png", plot = gg,units = "cm",height = 9, width = 15,dpi = 600, bg = "white")
# ggsave("figures_data_n_scripts/RankingRandomTFtargetsLists_terms.tiff", plot = gg,units = "cm",height = 9, width = 15,dpi = 600, bg = "white")


############ THE SAME BUT WITH FACET_WRAP

gg <- ggplot(data = topBiassedAnnotsDF, aes(y = RankingMedian)) +
  geom_pointrange(aes(x = as.numeric(size) + x_vals,
                      color = class, fill = class, 
                      ymin = RankingQ1, ymax = RankingQ3),
                  linetype='solid',
                  size = 0.001, fatten = 0.01) +
  geom_point(aes(x = as.numeric(size) + x_vals, y = RankingMedian,
                 color = class, fill = class),
             size = 0, stroke = 0) + 
  scale_x_discrete(labels = levels(topBiassedAnnotsDF$size))+
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 6)) +
  scale_fill_manual(values = color_pal,guide = 'none') +
  scale_color_manual(values = color_pal,guide = 'none') + 
  xlab("Size of TFs Lists") +
  ylab("Median Ranking Position") + 
  geom_text_repel(data = topBiassedAnnotsDF,
                  aes(x = size,
                      y = RankingMedian,
                      label = swr(term,15)),
                  seed = 343435,
                  min.segment.length = unit(0, "lines"),
                  segment.size = 0.2,
                  size = 2,
                  lineheight = 0.6,
                  box.padding = 0.5,
                  force_pull = 1) + 
  facet_wrap(~annotation,2,scales = "free")

ggsave("figures_data_n_scripts/RankingRandomTFtargetsLists_termsFACET.png", plot = gg,units = "cm",height = 9, width = 15,dpi = 600, bg = "white")
# ggsave("figures_data_n_scripts/RankingRandomTFtargetsLists_termsFACET.tiff", plot = gg,units = "cm",height = 9, width = 15,dpi = 600, bg = "white")


