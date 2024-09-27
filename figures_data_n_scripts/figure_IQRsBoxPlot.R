library(vroom)
library(ggplot2)
library(ggrepel)
#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(tidyverse)
library(colorspace)
library(ggrepel)
library(cowplot)
set.seed(1234)
library(stringr)
swr = function(string, nwrap=50) {paste(strwrap(string, width=nwrap), collapse="\n")}; swr = Vectorize(swr)

currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")
l_color <- lighten(colors_blind, 0.2)
d_color <- darken(colors_blind, 0.2)
color_pal <- c("Hypergeometric" = d_color, "Non-Central Hypergeometric" = l_color)

masterDF <- read.delim("random_lists_analysis/masterSimsResults.tsv") 
masterDF$annotation <- gsub('_', ' ', masterDF$annotation)
annotInfo <- read.delim("data/annotation_info_table.tsv")
annotInfo <- rbind(annotInfo,cbind(annotation_id="GO:0042493", term="response to xenobiotic stimulus"))

masterDF <- merge(annotInfo,masterDF,by.x ="annotation_id",by.y = "term",all.y = T)


# masterDF <- masterDF %>% pivot_longer(c(pTargetsF, pTargetsW), names_to = "method", values_to = "count") %>%
#   mutate(method = ifelse(method == "pTargetsF", "Hypergeometric", "Non-Central Hypergeometric"))

submasterDFQ1 <- masterDF %>% select(
  annotation_id, term, 
  Q1rTargetsF, Q1rTargetsW, 
  size,
  annotation
)

submasterDFQ3 <- masterDF %>% select(
  annotation_id,
  term,
  Q3rTargetsF, Q3rTargetsW,
  size,
  annotation
)

submasterDFRank <- masterDF %>% select(
  annotation_id, term, 
  rTargetsF, rTargetsW,
  size,
  annotation
)


submasterDFQ1 <- submasterDFQ1 %>% pivot_longer(c(Q1rTargetsF, Q1rTargetsW), names_to = "method", values_to = "RankingQ1") %>% mutate(method = ifelse(method == "Q1rTargetsF", "Hypergeometric", "Non-Central Hypergeometric"))
submasterDFQ3 <- submasterDFQ3 %>% pivot_longer(c(Q3rTargetsF, Q3rTargetsW), names_to = "method", values_to = "RankingQ3") %>% mutate(method = ifelse(method == "Q3rTargetsF", "Hypergeometric", "Non-Central Hypergeometric"))
submasterDFRank <- submasterDFRank %>% pivot_longer(c(rTargetsF, rTargetsW), names_to = "method", values_to = "RankingMedian") %>% mutate(method = ifelse(method == "rTargetsF", "Hypergeometric", "Non-Central Hypergeometric"))
#mutate(method = ifelse(method == "pTargetsF", "Fisher's Exact", "Wallenius")) ??? better this ???

final_masterDF <- merge(submasterDFRank,merge(submasterDFQ1,submasterDFQ3))



set.seed(98)
# masterDF <- final_masterDF %>% group_by(annotation) %>% mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>% ungroup() %>% mutate(db = factor(annotation, levels = names(colors_blind)), size = factor(size, levels=c(3,10,15,20)))


# Methods to differentiate
method_levels <- unique(final_masterDF$method)

topToShow <- 10
topBiassedAnnotsDF <- final_masterDF[final_masterDF$method == "Hypergeometric",] %>% arrange(RankingMedian) %>% group_by(annotation,size) %>% slice_head(n=topToShow)

selectedTerms <- paste0(topBiassedAnnotsDF$annotation_id,"_",topBiassedAnnotsDF$size)
topBiassedAnnotsDF <- final_masterDF[paste0(final_masterDF$annotation_id,"_",final_masterDF$size) %in% selectedTerms,] 

write.table(topBiassedAnnotsDF,file = "figures_data_n_scripts/RankingPlots_topBiassedAnnots.tsv",sep = "\t",quote = F,row.names = F,col.names = T)

topBiassedAnnotsDF$class <- paste0(topBiassedAnnotsDF$method,".",topBiassedAnnotsDF$annotation)
topBiassedAnnotsDF$colors <- color_pal[match(topBiassedAnnotsDF$class,names(color_pal))]

topBiassedAnnotsDF <- topBiassedAnnotsDF %>% mutate(annotation = factor(annotation, levels = names(colors_blind)), size=factor(size,levels = sort(unique(size)))) 

topBiassedAnnotsDF <- topBiassedAnnotsDF %>% group_by(method) %>% mutate(n = n(), x_vals = runif(n, -0.3, 0.3)) %>% ungroup() %>%  mutate(x_vals = ifelse(method == "Hypergeometric", x_vals -0.25, x_vals + 0.25))

annotRankingTops <- list()
for (ann in levels(topBiassedAnnotsDF$annotation)){
  sub_topBiassedAnnotsDF <- topBiassedAnnotsDF[topBiassedAnnotsDF$annotation == ann,]
  top_ann <- sub_topBiassedAnnotsDF %>% arrange(RankingMedian) %>% group_by(annotation,size,method) 
  annotRankingTops[[ann]] <- ggplot(data = sub_topBiassedAnnotsDF, aes(y = RankingMedian)) +
      geom_pointrange(aes(x = as.numeric(factor(size)) + x_vals,
                        color = class, fill = class, ymin = RankingQ1, ymax = RankingQ3), #,stroke=0),
                    linetype='solid', lwd=0.3, shape=18 ,size = 0.1, fatten = 2) +
    theme_classic() +
    theme(legend.position = "bottom", text = element_text(size = 6),
          panel.grid.major.y = element_line(color = "black",
                                            size = 0.1,
                                            linetype = 1)) +
    scale_fill_manual(values = color_pal,guide = 'none') +
    scale_color_manual(values = color_pal,guide = 'none') + 
    geom_text_repel(data = top_ann %>% slice_head(n=5),
      aes(x = as.numeric(factor(size)) + x_vals,
          y = RankingMedian,
          label = swr(term,15)),
          seed = 343435,
          min.segment.length = unit(0, "lines"),
          segment.size = 0.2,
          size = 1,
          lineheight = 0.6,
          box.padding = 0.5,
          force_pull = 1) + 
    xlab("Size of TFs Lists") +
    ylab("Median Ranking Position") + 
    scale_y_continuous(breaks =c(1, 25, 50, 100, round(max(sub_topBiassedAnnotsDF$RankingQ3))))
  #  + geom_boxplot(
  #   aes(
  #     x = size,
  #     color = method,
  #     color = after_scale(darken(color, .5, space = "HLS")),
  #     fill = after_scale(desaturate(lighten(color, .8), .4))
  #   ),
  #   outlier.shape = NA,
  #   alpha = 0.8,
  #   width = 0.6,
  #   linewidth = 0.08
  # )+
  # scale_x_discrete(labels = levels(masterDF_sub$size))+
}
gg <- plot_grid(plotlist = annotRankingTops)

ggsave("figures_data_n_scripts/RankingRandomTFtargetsLists_terms.png", plot = gg,units = "cm",height = 10, width = 16,dpi = 600, bg = "white")

ggsave("figures_data_n_scripts/RankingRandomTFtargetsLists_terms.tiff", plot = gg,units = "cm",height = 10, width = 16,dpi = 600, bg = "white",compression = "lzw")


############ THE SAME BUT WITH FACET_WRAP
# gg <- ggplot(data = topBiassedAnnotsDF, aes(y = RankingMedian)) +
#   geom_pointrange(aes(x = as.numeric(size) + x_vals,
#                       color = class, fill = class, 
#                       ymin = RankingQ1, ymax = RankingQ3),
#                   linetype='solid', shape=4,
#                   size = 1, fatten = .1) +
#   scale_x_discrete(labels = levels(topBiassedAnnotsDF$size))+
#   theme_classic() +
#   theme(legend.position = "bottom", text = element_text(size = 6)) +
#   scale_fill_manual(values = color_pal,guide = 'none') +
#   scale_color_manual(values = color_pal,guide = 'none') + 
#   xlab("Size of TFs Lists") +
#   ylab("Median Ranking Position") + 
#   geom_text_repel(data = topBiassedAnnotsDF %>% slice_head(n=3),
#                   aes(x = size, 
#                       y = RankingMedian,
#                       label = swr(term,15)),
#                   seed = 343435,
#                   min.segment.length = unit(0, "lines"),
#                   segment.size = 0.2,
#                   size = 0.6,
#                   lineheight = 0.6,
#                   box.padding = 0.5,
#                   force_pull = 1) + 
#   facet_wrap(~annotation,2,scales = "free")
# 
# ggsave("figures_data_n_scripts/RankingRandomTFtargetsLists_termsFACET.png", plot = gg,units = "cm",height = 9, width = 15,dpi = 600, bg = "white")
# ggsave("figures_data_n_scripts/RankingRandomTFtargetsLists_termsFACET.tiff", plot = gg,units = "cm",height = 9, width = 15,dpi = 600, bg = "white",compression = "lzw")
