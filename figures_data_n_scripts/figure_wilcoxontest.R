library(ggplot2)
library(ggrepel)
library(colorspace)
currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))

################################################################################
### 1 GENERATE DATA FOR PLOT
################################################################################
### Mann-Whitney U test (Wilcoxon rank-sum test)

annotationFiles <- c('data/GO_BP.tsv', 'data/KEGG.tsv', 'data/Reactome.tsv', 'data/WikiPathways.tsv'); annotationFile <- annotationFiles[1]
resultsDF <- c()
for (annotationFile in annotationFiles){
  annotation <- gsub('.tsv','',basename(annotationFile))
  annotationDF <- read.delim(annotationFile,header = T,sep = '\t')
  annotationDF <- annotationDF[annotationDF$organism == 9606,]
  collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')
  collectriTFsGRN <- collectriTFsGRN[!grepl("_",collectriTFsGRN$source),]
  
  collectriTFsGRN <- collectriTFsGRN[collectriTFsGRN$target_genesymbol %in% annotationDF$symbol, ]
  annotationDF <- annotationDF[annotationDF$symbol %in% collectriTFsGRN$target_genesymbol,]
  
  mergedData <- merge(annotationDF,collectriTFsGRN,by.x = "symbol", by.y="target_genesymbol")
  
  TFsperTargetAnnot <- as.data.frame(t(table(collectriTFsGRN$target_genesymbol)))[,2:3]
  colnames(TFsperTargetAnnot)[1] <- "target_genesymbol"
  
  for (annot in unique(annotationDF$annotation_id)){
    # annot <- unique(annotationDF$annotation_id)[1]
    TargetsAnnot <- unique(mergedData$symbol[mergedData$annotation_id == annot])
    annotTargetsPerTF <- TFsperTargetAnnot$Freq[TFsperTargetAnnot$target_genesymbol %in% TargetsAnnot]
    res <- wilcox.test(annotTargetsPerTF, TFsperTargetAnnot$Freq)
    resultsDF <- rbind(resultsDF,cbind(annot,res$p.value,annotation))
  }
}  
resultsDF <- as.data.frame(resultsDF)
colnames(resultsDF) <- c("annotation_id","wilcox.test","annotation")
write.table(resultsDF, file = 'data/observations/wilcox.test.tsv' ,sep = '\t',
            row.names = F, col.names = T)


################################################################################
### 2. PLOT DATA
################################################################################

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

ggsave("figures_data_n_scripts/figure_wilcoxon.tiff", plot = p,units = "cm",height = 6, width = 8,dpi = 500, compression = "lzw",bg = "white")

ggsave("figures_data_n_scripts/figure_wilcoxon.png", plot = p,units = "cm",height = 6, width = 8,dpi = 500, bg = "white")
