library(vroom)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(colorspace)
library(cowplot)

colors_blind <- c(
  "KEGG" = "#979A61",
  "GO BP" = "#3d91e0",
  "WikiPathways" = "#D974A0",
  "Reactome" = "#f1b620"
)

masterDF <- read.delim("data/EnrResults/masterSummary.tsv") 
masterDF <- masterDF %>% select(
  terms,
  TFsTargetsUni_HypergeomCount,
  size,
  annotation
)
masterDF$annotation <- gsub('_', ' ', masterDF$annotation)

set.seed(343435)
masterDF <- masterDF %>% group_by(annotation) %>% mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>%
  ungroup() %>% mutate(db = factor(annotation, levels = names(colors_blind)),
                       size = factor(size, levels=c(3,10,15,20)))

top_ann <- masterDF %>% arrange(desc(TFsTargetsUni_HypergeomCount)) %>% group_by(db,size) %>% slice_head(n=5)
top_ann <- top_ann %>% mutate(db = factor(db, levels = names(colors_blind))) %>%
  mutate(db_numeric = as.numeric(factor(db)))

gg <- ggplot(data = masterDF, aes(y = TFsTargetsUni_HypergeomCount)) +
  geom_point(aes(
    x = as.numeric(factor(size)) + x_vals,
    color = db,
    size = TFsTargetsUni_HypergeomCount
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
  facet_wrap(~db, scales = "free")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 6)) +
  scale_size(range = c(0.01, 2))+
  xlab("Size of Lists of TFs")+
  ylab("Mean(p-value < 0.05)")+
  scale_fill_manual(values = colors_blind) +
  scale_color_manual(values = colors_blind) +
  geom_text_repel(
    data = top_ann,
    aes(x = as.numeric(factor(size)) + x_vals, label = terms),
    seed = 343435,
    min.segment.length = unit(0, "lines"),
    segment.size = 0.2,
    size = 1.5
  )


ggsave("Final/figures_data_n_scripts/RandomTFsHypergeom.tiff", plot = gg,units = "cm",height = 12, width = 20,dpi = 300)


#### Targets

masterDF <- read.delim("data/EnrResults/masterSummary.tsv") %>% select(
  terms,
  Target_HypergeomCount,
  Target_NonCentralCount,
  size,
  annotation
)
masterDF$annotation <- gsub('_', ' ', masterDF$annotation)

masterDF <- masterDF %>% pivot_longer(c(Target_HypergeomCount, Target_NonCentralCount), names_to = "method", values_to = "count") %>%
  mutate(method = ifelse(method == "Target_HypergeomCount", "Hypergeometrics", "NonCentral Hypergeometrics"))

set.seed(343435)
masterDF <- masterDF %>% group_by(annotation) %>% mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>%
  ungroup() %>% mutate(db = factor(annotation, levels = names(colors_blind)),
                       size = factor(size, levels=c(3,10,15,20)))

# Create a new column combining db and method for coloring
masterDF <- masterDF %>%
  mutate(db_method = paste(db, method, sep = "_"))

# Methods to differentiate
method_levels <- unique(masterDF$method)
ann <- "KEGG"
splitted_plots <- lapply(levels(masterDF$db), function(ann){
  masterDF_sub <- masterDF %>% filter(db == ann)
  color_pal <- colors_blind[[ann]]
  l_color <- lighten(color_pal, 0.2)
  d_color <- darken(color_pal, 0.2)
  color_pal <- c("Hypergeometrics" = l_color, "NonCentral Hypergeometrics" = d_color)
  
  masterDF_sub <- masterDF_sub %>% group_by(method) %>% 
    mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>% ungroup() %>%
    mutate(x_vals = ifelse(method == "Hypergeometrics", x_vals -0.15, x_vals + 0.15))
  
  ggplot(data = masterDF_sub, aes(y = count)) +
    geom_point(aes(
      x = as.numeric(factor(size)) + x_vals,
      color = method,
      size = count
    )) +
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
    scale_x_discrete(labels = levels(masterDF$size))+
    facet_wrap(~db, scales = "free")+
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
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal)+
    scale_size(range = c(0.01, 2), guide = "none")
})

gg <- plot_grid(plotlist = splitted_plots)

gg <- grid.arrange(splitted_plots, left="JAMON")

ggsave("Final/figures_data_n_scripts/RandomTFtargetsLists.tiff", plot = gg,units = "cm",height = 6, 
       width = 8,dpi = 600, bg = "white")
