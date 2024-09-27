############# Plots Collectri Targets and TFs distribution in Annotations DBs 
"
This plot answers the questions:
1- How is the distribution of targets per annotation? 
2- How is the distribution of TFs per annotation?
"
currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))

library(vroom)
library(ggplot2)
library(ggrepel)
#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(ggpp)
library(tidyverse)
library(ggdist)
library(colorspace)
library(ggrepel)

################################################################################
### 1 GENERATE DATA FOR PLOT
################################################################################
### CollecTri TF - Target Distribution - Per Annotation DB

collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv'); annotationFile <- annotationFiles[1]
TFsPerTargetdf_full <- c()
for (annotationFile in annotationFiles){
  annotation <- gsub('.tsv','',basename(annotationFile))
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == 9606,]
  collectriTFsGRN_filtered <- collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% annotationDF$symbol,]
  TFsPerTargetdf <- as.data.frame(table(collectriTFsGRN_filtered$target_genesymbol))
  colnames(TFsPerTargetdf) <- c('target','numberOfTFsAssociated')
  TFsPerTargetdf$db <- annotation
  TFsPerTargetdf_full <- rbind(TFsPerTargetdf_full, TFsPerTargetdf)
}
View(TFsPerTargetdf_full)

write.table(TFsPerTargetdf_full, file = 'data/observations/TFsPerTargetdf.tsv',
            sep = '\t', row.names = F, col.names = T)

################################################################################
### 2. PLOT DATA
################################################################################

TFsNtargets <- read.delim("data/observations/TFsPerTargetdf.tsv")
TFsNtargets$db <- gsub('_',' ',TFsNtargets$db)

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

bx_limits <- TFsNtargets %>%
  group_by(db) %>%
  summarise(stats = list(boxplot_limits(numberOfTFsAssociated))) %>%
  unnest_wider(stats)

n_exceed_limits <- TFsNtargets %>% 
  group_by(db) %>% 
  mutate(total = n()) %>%
  ungroup() %>%
  inner_join(bx_limits) %>%
  filter(numberOfTFsAssociated > upper_whisker) %>% 
  group_by(db) %>%
  mutate(filtered = n()) %>%
  ungroup() %>%
  select(db, total, filtered, upper_whisker) %>%
  unique() %>%
  mutate(freq = filtered / total)

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

point_coords <- c()

set.seed(1234)
TFsNtargets <- TFsNtargets %>% group_by(db) %>% mutate(n =n(), y_vals = runif(n, -0.2, 0.2)) %>% 
  mutate(db = factor(db, levels = names(colors_blind))) %>% ungroup() %>%
  mutate(db_numeric = as.numeric(factor(db)))

n_exceed_limits <- n_exceed_limits %>% mutate(db = factor(db, levels = names(colors_blind))) %>%
  mutate(db_numeric = as.numeric(factor(db)))

top_genes <- TFsNtargets %>% arrange(desc(numberOfTFsAssociated)) %>% group_by(db) %>% slice_head(n=5)
top_genes <- top_genes %>% mutate(db = factor(db, levels = names(colors_blind))) %>%
  mutate(db_numeric = as.numeric(factor(db)))


gg <- ggplot(data = TFsNtargets, aes(x = numberOfTFsAssociated)) +
  geom_point(aes(
    y = as.numeric(factor(db)) + y_vals,
    color = db,
    size = numberOfTFsAssociated
  )) +
  geom_boxplot(
    aes(
      y = db,
      color = db,
      color = after_scale(darken(color, .5, space = "HLS")),
      fill = after_scale(desaturate(lighten(color, .8), .4))
    ),
    outlier.shape = NA,
    alpha = 0.8,
    width = 0.5,
    linewidth = 0.5
  ) +
  scale_y_discrete(labels = levels(TFsNtargets$db)) +
  labs(x = "Number of TFs Associated to each Target Gene") +
  theme_classic() +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  scale_size(range = c(0.01, 2)) +
  stat_boxplot(
    geom = 'errorbar',
    aes(
      y = db,
      color = db,
      color = after_scale(darken(color, .5, space = "HLS")),
      
    ),
    width = .3,
    linewidth = 0.5
  ) +
  geom_segment(
    data = n_exceed_limits,
    aes(
      y = as.numeric(factor(db)) + 0.3,
      x = upper_whisker,
      yend = as.numeric(factor(db)) + 0.3,
      xend = upper_whisker + 20
    ),
    arrow = arrow(length = unit(0.1, "cm"))
  ) +
  geom_text(
    data = n_exceed_limits,
    aes(
      y = as.numeric(factor(db)) + 0.3,
      label = paste0("Outlier Genes: ", filtered, " (", 100 * round(freq, 4), "%)"),
      x = upper_whisker + 25,
      hjust = 0
    ),
    size = 2.5
  ) +
  scale_fill_manual(values = colors_blind) +
  scale_color_manual(values = colors_blind) +
  geom_text_repel(
    data = top_genes,
    aes(y = as.numeric(factor(db)) + y_vals, label = target),
    seed = 1234,
    min.segment.length = unit(0, "lines"),
    segment.size = 0.2,
    size = 2
  )

ggsave("figures_data_n_scripts/collectriTargetsNtfsPerAnnDB.tiff", plot = gg,units = "cm",height = 12, width = 20,dpi = 500, compression = "lzw",bg = "white")

ggsave("figures_data_n_scripts/collectriTargetsNtfsPerAnnDB.png", plot = gg,units = "cm",height = 12, width = 20,dpi = 500, bg = "white")
