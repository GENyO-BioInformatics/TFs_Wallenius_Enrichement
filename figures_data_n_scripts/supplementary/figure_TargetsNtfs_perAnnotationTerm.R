library(vroom)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(colorspace)
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

TFsNtargetsAnnotaDist <- vroom("data/observations/TFsNtargetsPerAnnotationDistribution.tsv")
TFsNtargetsAnnotaDist$annotationDB <- gsub('_',' ',TFsNtargetsAnnotaDist$annotationDB)

TFsNtargetsAnnotaDist$annotationDB <- factor(TFsNtargetsAnnotaDist$annotationDB, levels = names(colors_blind))

TFsNtargetsAnnotaDist <- TFsNtargetsAnnotaDist %>% group_by(annotationDB) %>% mutate(ranking = rank(-relativeRisk))

ggplot(TFsNtargetsAnnotaDist, aes(y = relativeRisk, x = annotationDB))+
  geom_jitter(aes(fill = annotationDB, color = after_scale(darken(fill, .2, space = "HLS")),
                 fill = after_scale(lighten(fill, .2, space = "HLS"))), alpha = 0.5, shape = 21,
              width = 0.25)+
  geom_boxplot(aes(color = annotationDB,
                   color = after_scale(darken(color, .5, space = "HLS")),
                   fill = after_scale(desaturate(lighten(color, .8), .4))
                   ),
    outlier.shape = NA,
    alpha = 0.8,
    width = 0.5,
    linewidth = 0.5)+
  scale_fill_manual(values = colors_blind)+
  scale_color_manual(values = colors_blind)+
  labs(y = "Relative Risk")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_log10()

ggplot(TFsNtargetsAnnotaDist, aes(x = relativeRisk))+
  geom_density()+
  facet_wrap(~annotationDB, scales = "free")

ggsave("Figures/TargetGenes_vs_AverageTFs_perAnnotationTerm.tiff", plot = gg,units = "cm",height = 12, width = 20,dpi = 300)
