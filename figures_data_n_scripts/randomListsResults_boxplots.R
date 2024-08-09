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

masterDF <- read.delim("data/EnrResults/masterSummary.tsv") %>% select(
  terms,
  TFs_HypergeomCount,
  size,
  annotation
)
masterDF$annotation <- gsub('_', ' ', masterDF$annotation)

set.seed(343435)
masterDF <- masterDF %>% group_by(annotation) %>% mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>%
  ungroup() %>% mutate(db = factor(annotation, levels = names(colors_blind)),
                       size = factor(size, levels=c(3,10,15,20)))

top_ann <- masterDF %>% arrange(desc(TFs_HypergeomCount)) %>% group_by(db,size) %>% slice_head(n=5)
top_ann <- top_ann %>% mutate(db = factor(db, levels = names(colors_blind))) %>%
  mutate(db_numeric = as.numeric(factor(db)))

gg <- ggplot(data = masterDF, aes(y = TFs_HypergeomCount)) +
  geom_point(aes(
    x = as.numeric(factor(size)) + x_vals,
    color = db,
    size = TFs_HypergeomCount
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

ggsave("TFsSEA_boxplots.tiff", plot = gg,units = "cm",height = 12, width = 20,dpi = 300)

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
  mutate(method = ifelse(method == "Target_HypergeomCount", "Hypergeometric", "Non-Central Hypergeometric"))

set.seed(343435)
masterDF <- masterDF %>% group_by(annotation) %>% mutate(n = n(), x_vals = runif(n, -0.1, 0.1)) %>%
  ungroup() %>% mutate(db = factor(annotation, levels = names(colors_blind)),
                       size = factor(size, levels=c(3,10,15,20)))

# Create a new column combining db and method for coloring
masterDF <- masterDF %>%
  mutate(db_method = paste(db, method, sep = "_"))

# Generate the color palette using colorspace functions
generate_colors <- function(base_colors, method_levels) {
  pal <- setNames(rep("", length(base_colors) * length(method_levels)), 
                  paste(rep(names(base_colors), each = length(method_levels)), 
                        method_levels, sep = "_"))
  
  for (db in names(base_colors)) {
    base_color <- base_colors[db]
    for (method in method_levels) {
      color_name <- paste(db, method, sep = "_")
      if (method == "Hypergeometric") {
        pal[color_name] <- darken(base_color, amount = 0.3)  # Darker shade
      } else {
        pal[color_name] <- lighten(base_color, amount = 0.3)  # Lighter shade
      }
    }
  }
  pal
}

# Methods to differentiate
method_levels <- unique(masterDF$method)

splitted_plots <- lapply(levels(masterDF$db), function(ann){
  masterDF_sub <- masterDF %>% filter(db == ann)
  color_pal <- colors_blind[[ann]]
  l_color <- lighten(color_pal, 0.2)
  d_color <- darken(color_pal, 0.2)
  color_pal <- c("Hypergeometric" = d_color, "Non-Central Hypergeometric" = l_color)
  ggplot(data = masterDF_sub, aes(x = size, y = count)) +
    geom_point(aes(
      color = method,
      color = after_scale(darken(color, .5, space = "HLS")),
      fill = after_scale(desaturate(lighten(color, .8), .4)),
      size = count),
    position = position_jitterdodge()) +
    scale_size(range = c(0.001, 0.1)) +
    geom_boxplot(
      aes(color = method,
        color = after_scale(darken(color, .5, space = "HLS")),
        fill = after_scale(desaturate(lighten(color, .8), .4))
      ),
      outlier.size = 0.05,
      outlier.shape = NA,
      alpha = 0.8,
      width = 0.6,
      linewidth = 0.08) +
    scale_x_discrete(labels = levels(masterDF$size))+
    theme_classic()+
    theme(legend.position = "bottom",
          legend.key.size = unit(0.1,"cm"),
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
    scale_color_manual(values = color_pal)
})


gg <- plot_grid(plotlist = splitted_plots)

ggsave("TFsTargetsSEA_boxplots.tiff", plot = gg,units = "cm",height = 6, width = 8,dpi = 500,
       compression = "lzw",bg = "white")
