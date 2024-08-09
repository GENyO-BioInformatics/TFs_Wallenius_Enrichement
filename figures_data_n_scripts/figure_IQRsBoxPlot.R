############# Plots Collectri Targets and TFs distribution in Annotations DBs 
"
This plot answers the questions:
1- How is the distribution of targets per annotation? 
2- How is the distribution of TFs per annotation?

"

library(vroom)
library(ggplot2)
library(ggrepel)
#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(ggpp)
library(tidyverse)
library(ggdist)
library(colorspace)
library(ggrepel)

TargetsSEA <- read.delim("data/EnrResults/TargetsSEAresultsMelted.tsv")
TargetsSEA$annotation <- gsub('_',' ',TargetsSEA$annotation)

boxplot_limits <- function(data) {
  box_stats <- boxplot.stats(data)$stats
  list(
    lower_whisker = box_stats[1],
    q1 = box_stats[2],
    median = box_stats[3],
    q3 = box_stats[4],
    upper_whisker = box_stats[5]
  )
}

bx_limits <- TargetsSEA %>%
  group_by(annotation) %>%
  summarise(stats = list(boxplot_limits(RankingIQR))) %>%
  unnest_wider(stats)

n_exceed_limits <- TargetsSEA %>% 
  group_by(annotation) %>% 
  mutate(total = n()) %>%
  ungroup() %>%
  inner_join(bx_limits) %>%
  filter(RankingIQR > upper_whisker) %>% 
  group_by(annotation) %>%
  mutate(filtered = n()) %>%
  ungroup() %>%
  select(annotation, total, filtered, upper_whisker) %>%
  unique() %>%
  mutate(freq = filtered / total)

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

point_coords <- c()

set.seed(1234)
TargetsSEA <- TargetsSEA %>% group_by(annotation) %>% mutate(n =n(), y_vals = runif(n, -0.2, 0.2)) %>% 
  mutate(annotation = factor(annotation, levels = names(colors_blind))) %>% ungroup() %>%
  mutate(db_numeric = as.numeric(factor(annotation)))

n_exceed_limits <- n_exceed_limits %>% mutate(annotation = factor(annotation, levels = names(colors_blind))) %>%
  mutate(db_numeric = as.numeric(factor(annotation)))

top_genes <- TargetsSEA %>% arrange(desc(RankingIQR)) %>% group_by(annotation) %>% slice_head(n=5)
top_genes <- top_genes %>% mutate(annotation = factor(annotation, levels = names(colors_blind))) %>%
  mutate(db_numeric = as.numeric(factor(annotation)))


gg <- ggplot(data = TargetsSEA, aes(x = RankingIQR)) +
  geom_point(aes(
    y = as.numeric(factor(annotation)) + y_vals,
    color = annotation,
    size = RankingIQR
  )) +
  geom_boxplot(
    aes(
      y = annotation,
      color = annotation,
      color = after_scale(darken(color, .5, space = "HLS")),
      fill = after_scale(desaturate(lighten(color, .8), .4))
    ),
    outlier.shape = NA,
    alpha = 0.8,
    width = 0.5,
    linewidth = 0.5
  ) +
  scale_y_discrete(labels = levels(TargetsSEA$annotation)) +
  labs(x = "IQR Associated to Each Annotation") +
  theme_classic() +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  scale_size(range = c(0.01, 2)) +
  stat_boxplot(
    geom = 'errorbar',
    aes(
      y = annotation,
      color = annotation,
      color = after_scale(darken(color, .5, space = "HLS")),
    ),
    width = .3,
    linewidth = 0.5
  ) +
  geom_segment(
    data = n_exceed_limits,
    aes(
      y = as.numeric(factor(annotation)) + 0.3,
      x = upper_whisker,
      yend = as.numeric(factor(annotation)) + 0.3,
      xend = upper_whisker + 20
    ),
    arrow = arrow(length = unit(0.1, "cm"))
  ) +
  scale_fill_manual(values = colors_blind) +
  scale_color_manual(values = colors_blind) +
  geom_text_repel(
    data = top_genes,
    aes(y = as.numeric(factor(annotation)) + y_vals, label = annotation_id),
    seed = 1234,
    min.segment.length = unit(0, "lines"),
    segment.size = 0.2,
    size = 2
  )

ggsave("TargetsSEAresults_IQRBoxplot.tiff", plot = gg,units = "cm",height = 12, width = 20,dpi = 500,
       compression = "lzw",bg = "white")

