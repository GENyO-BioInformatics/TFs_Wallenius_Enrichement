####################### Plots Annotations Coverage By Collectri Evidence 
"
This plot answers two questions:
1- Which % of each annotation database (annotations / total annotation * 100) 
   is covered in function of the evidence level in Collectri target genes?
2- Which % of each annotation database (annotations / total annotation * 100) 
   is covered in function of the evidence level in Collectri TFs?
"

library(vroom)
library(ggplot2)

annotationsCoveragePerEvidence <- vroom("data/observations/annotationsCoveragePerEvidence.tsv")
annotationsCoveragePerEvidence$annotation <- gsub('_',' ',annotationsCoveragePerEvidence$annotation)
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

# byN_references NO APORTA NADA
# byN_references <- annotationsCoveragePerEvidence[annotationsCoveragePerEvidence$evidence == 'n_references',]
byCuration_effortDF <- annotationsCoveragePerEvidence[annotationsCoveragePerEvidence$evidence == 'curation_effort',]

newDF_targets <- byCuration_effortDF
newDF_targets$propor_TFsAnnotationsInEvidence <- 'Targets'

newDF_TFs <- byCuration_effortDF
newDF_TFs$propor_TargetsAnnotationsInEvidence <- newDF_TFs$propor_TFsAnnotationsInEvidence
newDF_TFs$propor_TFsAnnotationsInEvidence <- 'TFs'

byCuration_effortDF <- rbind(newDF_targets,newDF_TFs)
colnames(byCuration_effortDF)[colnames(byCuration_effortDF) == 'propor_TFsAnnotationsInEvidence'] <- 'Entity'
colnames(byCuration_effortDF)[colnames(byCuration_effortDF) == 'propor_TargetsAnnotationsInEvidence'] <- 'AnnotationDBCoverage'

coveragePlot <- ggplot(byCuration_effortDF, 
       aes(x=minEvidenceScore, y=AnnotationDBCoverage,
           color=annotation, linetype=Entity)) +
  scale_color_manual(values = colors_blind, breaks = names(colors_blind), name = 'Annotation\nDatabase') +
  geom_line(size=0.3) +
  theme(panel.background = element_blank(),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        legend.key.size = unit(0.2,"cm"),
        legend.key = element_rect(fill = "white"),
        plot.margin = margin(l = 0, unit = "cm"))+
  ylab("Annotation coverage (%)")+
  xlab("Evidence Level") +
  labs(fill = "Annotation")

coveragePlot

ggsave("figure_collectriAnnotationsCoverage.tiff", plot = coveragePlot,units = "cm",height = 6, width = 8,dpi = 500,
       compression = "lzw",bg = "white")

