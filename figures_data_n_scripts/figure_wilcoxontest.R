library(ggplot2)
library(ggrepel)
library(colorspace)


colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "Reactome" = "#f1b620",
                  "WikiPathways" = "#D974A0")

wilcoxTestDF <- read.delim("data/observations/wilcox.test.tsv")
wilcoxTestDF$wilcox.test
wilcoxTestDF$annotation_id
wilcoxTestDF$annotation <- gsub("_"," ",wilcoxTestDF$annotation)
wilcoxTestDF$annotation <- factor(wilcoxTestDF$annotation,levels = names(colors_blind))

brks <- seq(min(wilcoxTestDF$wilcox.test), max(wilcoxTestDF$wilcox.test), length.out = 30)

p <- ggplot()+
  geom_histogram(data = wilcoxTestDF, aes(x = wilcox.test, fill = annotation), breaks = brks, color = "gray20", linewidth = 0.2)+
  theme_classic(base_size = 4.5,base_line_size = 0.2)+
  theme(legend.position = "none",
        strip.background = element_blank())+
  ylab("Frequency")+
  xlab("Mann-Whitney U test p-value")+
  scale_fill_manual(values = colors_blind)+
  facet_wrap(~annotation,scales = "free_y")

ggsave("figure_wilcoxon.tiff", plot = p,units = "cm",height = 6, width = 8,dpi = 500,
       compression = "lzw",bg = "white")

ggsave("figure_wilcoxon.png", plot = p,units = "cm",height = 6, width = 8,dpi = 500,
       bg = "white")
