####################### Plots Annotations Coverage By Collectri Evidence 
"
This plot answers two questions:
1- Which % of each annotation database (annotations / total annotation * 100) 
   is covered in function of the evidence level in Collectri target genes?
2- Which % of each annotation database (annotations / total annotation * 100) 
   is covered in function of the evidence level in Collectri TFs?
"

currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))
library(vroom)
library(ggplot2)

################################################################################
### 1 GENERATE DATA FOR PLOT
################################################################################

collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv'); annotationFile <- annotationFiles[1]

evidences <- c('n_references','curation_effort'); evidence <- evidences[1]
annotationsCoveragePerEvidence <- c()
for (evidence in evidences){
  for (annotationFile in annotationFiles){
    annotation <- gsub('.tsv','',basename(annotationFile))
    outfile <- gsub('.tsv',paste0('_CollectriCoverage_',evidence,'.tsv'),annotationFile)
    outfile <- gsub('data','data/observations',outfile)
    
    annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
    annotationDF <- annotationDF[annotationDF$organism == 9606,]
    
    totalAnnotations <- length(unique(annotationDF$annotation_id))
    for (evidenceScore in sort(unique(collectriTFsGRN[,evidence]))){
      selection <- which(collectriTFsGRN[,evidence] >= evidenceScore)
      TFsAccepted <- unique(collectriTFsGRN$source_genesymbol[selection])
      TargetsAccepted <- unique(collectriTFsGRN$target_genesymbol[selection])
      nTFsAccepted <- length(TFsAccepted)
      nTargetsAccepted <- length(TargetsAccepted)
      TFsAnnotationsInEvidence <- length(unique(annotationDF$annotation_id[annotationDF$symbol %in% TFsAccepted]))
      TargetsAnnotationsInEvidence <- length(unique(annotationDF$annotation_id[annotationDF$symbol %in% TargetsAccepted]))
      propor_TFsAnnotationsInEvidence <- TFsAnnotationsInEvidence / totalAnnotations * 100
      propor_TargetsAnnotationsInEvidence <- TargetsAnnotationsInEvidence / totalAnnotations * 100
      
      annotationsCoveragePerEvidence <- rbind(annotationsCoveragePerEvidence,
                                              cbind(annotation=annotation,
                                                    evidence = evidence,
                                                    minEvidenceScore = evidenceScore,
                                                    nTFsAccepted = nTFsAccepted, 
                                                    TFsAnnotationsInEvidence = TFsAnnotationsInEvidence,
                                                    propor_TFsAnnotationsInEvidence = propor_TFsAnnotationsInEvidence,
                                                    nTargetsAccepted = nTargetsAccepted,
                                                    TargetsAnnotationsInEvidence = TargetsAnnotationsInEvidence,
                                                    propor_TargetsAnnotationsInEvidence= propor_TargetsAnnotationsInEvidence))
    }
  }
}
write.table(annotationsCoveragePerEvidence, file = 'data/observations/annotationsCoveragePerEvidence.tsv' ,sep = '\t',
            row.names = F, col.names = T)
View(annotationsCoveragePerEvidence)

################################################################################
### 2. PLOT DATA
################################################################################

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
  xlab("Curation Effort") +
  labs(fill = "Annotation")

ggsave("figures_data_n_scripts/figure_collectriAnnotationsCoverage.tiff", plot = coveragePlot,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")

ggsave("figures_data_n_scripts/figure_collectriAnnotationsCoverage.png", plot = coveragePlot,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")

