library(vroom)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(colorspace)
library(cowplot)

currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))

colors_blind <- c(
  "KEGG" = "#979A61",
  "GO BP" = "#3d91e0",
  "WikiPathways" = "#D974A0",
  "Reactome" = "#f1b620"
)

library(stringr)
swr = function(string, nwrap=50) {paste(strwrap(string, width=nwrap), collapse="\n")}; swr = Vectorize(swr)

masterDF <- read.delim("data/EnrResults/masterSummary.tsv") 
masterDF$annotation <- gsub('_', ' ', masterDF$annotation)
annotInfo <- read.delim("data/annotation_info_table.tsv")
masterDF <- merge(annotInfo,masterDF,by.x ="annotation_id",by.y = "terms",all.y = T)


# BOX PLOTS VISUALIZATIONS 

# TFs SEA ########################
masterDF <- masterDF %>% select(
  annotation_id,
  term,
  TFsTargetsUni_HypergeomCount,
  size,
  annotation
)

set.seed(343435) # Used by gg_text_repel @RAUL confirm this is why you set a seed

masterDF <- masterDF %>% group_by(annotation) %>% mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>%
  ungroup() %>% mutate(db = factor(annotation, levels = names(colors_blind)),
                       size = factor(size, levels=c(3,10,15,20)))

gg <- ggplot(data = masterDF, aes(y = TFsTargetsUni_HypergeomCount)) +
  geom_point(aes(
    x = as.numeric(factor(size)) + x_vals,
    color = db,
    fill = db,
    size = TFsTargetsUni_HypergeomCount,
    alpha = 0.6,
    stroke = 0
  )) +
  geom_boxplot(
    aes(
      x = size,
      color = db,
      color = after_scale(darken(color, .5, space = "HLS")),
      fill = after_scale(desaturate(lighten(color, .8), .4))
    ),
    outlier.shape = NA,
    alpha = 0.8,
    width = 0.3,
    linewidth = 0.5
  ) +
  scale_x_discrete(labels = levels(masterDF$size))+
  facet_wrap(~db, scales = "fixed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 6)) +
  scale_size(range = c(0.01, 2))+
  xlab("Size of TFs Lists")+
  ylab("Frequency(p-value < 0.05)")+
  scale_fill_manual(values = colors_blind) +
  scale_color_manual(values = colors_blind) +
  scale_alpha(guide = 'none')

top_ann <- masterDF %>% arrange(desc(TFsTargetsUni_HypergeomCount)) %>% group_by(db,size) 
top_ann <- top_ann %>% mutate(db = factor(db, levels = names(colors_blind))) %>%  mutate(db_numeric = as.numeric(factor(db)))

topToShow <- 3
# gg_ids <- gg + geom_text_repel(
#     data = top_ann %>% slice_head(n=topToShow),
#     aes(x = as.numeric(factor(size)) + x_vals, label = annotation_id),
#     seed = 343435,
#     min.segment.length = unit(0, "lines"),
#     segment.size = 0.2,
#     size = 2
# )
# 
# ggsave("figures_data_n_scripts/RandomTFsHypergeom_ids.tiff", plot = gg_ids,units = "cm",height = 12, width = 20,dpi = 300)
# 
# ggsave("figures_data_n_scripts/RandomTFsHypergeom_ids.png", plot = gg_ids,units = "cm",height = 12, width = 20,dpi = 300)

gg_terms <- gg + geom_text_repel(
  data = top_ann %>% slice_head(n=topToShow),
  aes(x = as.numeric(factor(size)) + x_vals, label = swr(term,15)),
  seed = 343435,
  min.segment.length = unit(0, "lines"),
  segment.size = 0.2,
  size = 1.5,
  lineheight = 0.6,
  box.padding = 0.2,
  force_pull = 2
)

ggsave("figures_data_n_scripts/RandomTFsHypergeom_terms.tiff", plot = gg_terms,units = "cm",height = 10, width = 15,dpi = 300)

ggsave("figures_data_n_scripts/RandomTFsHypergeom_terms.png", plot = gg_terms,units = "cm",height = 10, width = 15,dpi = 300)



# TFs Targets SEA ##############

masterDF <- read.delim("data/EnrResults/masterSummary.tsv") 
masterDF$annotation <- gsub('_', ' ', masterDF$annotation)
annotInfo <- read.delim("data/annotation_info_table.tsv")
masterDF <- merge(annotInfo,masterDF,by.x ="annotation_id",by.y = "terms",all.y = T)

masterDF <- masterDF %>% select(
  annotation_id,
  term,
  Target_HypergeomCount,
  Target_NonCentralCount,
  size,
  annotation
)

masterDF <- masterDF %>% pivot_longer(c(Target_HypergeomCount, Target_NonCentralCount), names_to = "method", values_to = "count") %>% 
  mutate(method = ifelse(method == "Target_HypergeomCount", "Hypergeometric", "Non-Central Hypergeometric"))
  #mutate(method = ifelse(method == "Target_HypergeomCount", "Fisher's Exact", "Wallenius")) ??? better this ???

set.seed(343435)
masterDF <- masterDF %>% group_by(annotation) %>% mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>% ungroup() %>% mutate(db = factor(annotation, levels = names(colors_blind)),
                       size = factor(size, levels=c(3,10,15,20)))

# Create a new column combining db and method for coloring
masterDF <- masterDF %>%  mutate(db_method = paste(db, method, sep = "_"))

# Methods to differentiate
method_levels <- unique(masterDF$method)

ann <- "GO BP"
topToShow <- 3
TFstargetsIds <- list()
TFstargetsTerms <- list()
for (ann in levels(masterDF$db)){
  masterDF_sub <- masterDF %>% filter(db == ann)
  color_pal <- colors_blind[[ann]]
  l_color <- lighten(color_pal, 0.2)
  d_color <- darken(color_pal, 0.2)
  color_pal <- c("Hypergeometric" = l_color, "Non-Central Hypergeometric" = d_color)
  
  masterDF_sub <- masterDF_sub %>% group_by(method) %>% 
    mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>% ungroup() %>%
    mutate(x_vals = ifelse(method == "Hypergeometric", x_vals - 0.15, x_vals + 0.15))
  
  gg <- ggplot(data = masterDF_sub, aes(y = count)) +
    geom_point(aes(
      x = as.numeric(factor(size)) + x_vals,
      color = method,
      fill = method,
      size = count,
      alpha = 0.6,
      stroke = 0)) +
    scale_size(range = c(0.01, 1), guide = "none") + 
    geom_boxplot(
      aes(
        x = size,
        color = method,
        color = after_scale(darken(color, .5, space = "HLS")),
        fill = after_scale(desaturate(lighten(color, .8), .4))
      ),
      outlier.shape = NA,
      alpha = 0.8,
      width = 0.6,
      linewidth = 0.08
    )+
    scale_x_discrete(labels = levels(masterDF_sub$size))+
    theme_classic()+
    theme(legend.position = "bottom",
          legend.key.size = unit(0.01,"cm"),
          legend.title = element_blank(),
          axis.title = element_blank(),
          text = element_text(size = 3),
          axis.line = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1),
          legend.background = element_blank(),
          legend.text = element_text(size = 2.5),
          legend.box.spacing = unit(0, "pt"),
          plot.margin = margin(1, 2, -1.5, 2),
          strip.background = element_blank())+
    scale_alpha(guide = 'none')+
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal)
  
  top_ann <- masterDF_sub %>% arrange(desc(count)) %>% group_by(db,size,method) 
  top_ann <- top_ann %>% mutate(db = factor(db, levels = names(colors_blind))) %>%  mutate(db_numeric = as.numeric(factor(db)))
  
  # gg_ids <- gg + geom_text_repel(
  #   data = top_ann %>% slice_head(n=topToShow),
  #   aes(x = as.numeric(factor(size)) + x_vals, label = annotation_id),
  #   seed = 343435,
  #   min.segment.length = unit(0, "lines"),
  #   segment.size = 0.2,
  #   size = .6
  # )
  
  gg_terms <- gg + geom_text_repel(
    data = top_ann %>% slice_head(n=topToShow),
    aes(x = as.numeric(factor(size)) + x_vals, label = swr(term,15)),
    seed = 8,
    min.segment.length = unit(0, "lines"),
    segment.size = 0.2,
    size = 1,
    lineheight = 0.65,
    box.padding = 0.2)
  
  TFstargetsIds[[ann]] <- gg_ids
  TFstargetsTerms[[ann]] <- gg_terms
}
  
gg <- plot_grid(plotlist = TFstargetsTerms)
ggsave("figures_data_n_scripts/RandomTFtargetsLists_terms.png", plot = gg,units = "cm",height = 9, width = 15,dpi = 600, bg = "white")
ggsave("figures_data_n_scripts/RandomTFtargetsLists_terms.tiff", plot = gg,units = "cm",height = 9, width = 15,dpi = 600, bg = "white")

# gg <- plot_grid(plotlist = TFstargetsIds)
# ggsave("figures_data_n_scripts/RandomTFtargetsLists_ids.png", plot = gg,units = "cm",height = 9, width = 8,dpi = 600, bg = "white")
# ggsave("figures_data_n_scripts/RandomTFtargetsLists_ids.tiff", plot = gg,units = "cm",height = 9, width = 8,dpi = 600, bg = "white")


