#
# Currently ranking is calculated based on raw pvalue.
# Should we use the adjusted pvalue? 
# When donw the rank function to produce extreme rankings and the plots result weird
#

library(vroom)
library(glue)
library(tidyverse)
library(parallel)
library(ggplot2)
library(sdamr)
library(ggrepel)
library(cowplot)
library(stringr)
library(data.table)
library(colorspace)
library(ggpubr)
library(gtable)
library(lemon)
library(stringr)
swr = function(string, nwrap=50) {paste(strwrap(string, width=nwrap), collapse="\n")}; swr = Vectorize(swr)

currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")
annotations <- names(colors_blind)


list_files <- list.files("use_case/results",pattern = "*_targetsF.tsv",full.names = T)
fileName <- list_files[1]
fullContent <- rbindlist(lapply(list_files, function(fileName){
  fileNameCheck <- gsub("GO_BP","GO BP",fileName)
  basenameFile <- basename(fileNameCheck)
  use_case <- gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[1])
  annotation <- gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[2])
  content <- vroom(fileName)
  names(content)[2] <- "pvalue"
  contentF <- content %>% mutate(method = "Hypergeometric",
                                pvalueAdj = p.adjust(pvalue),
                                ranking = rank(pvalueAdj)) # ranking = rank(pvalueAdj)) 
  contentF$ranking[is.na(contentF$pvalue)] <- NA
  
  content <- vroom(gsub("targetsF","targetsW",fileName))
  names(content)[2] <- "pvalue"
  contentW <- content %>% mutate(method = "Non-Central Hypergeometric",
                                pvalueAdj = p.adjust(pvalue), #
                                ranking = rank(pvalue)) #ranking = rank(pvalueAdj))
  contentW$ranking[is.na(contentW$pvalue)] <- NA
  
  mergedDF <- merge(contentF,contentW,suffixes = c("_F","_W"),by = "term")
  mergedDF <- mergedDF %>% mutate(method = "Non-Central Hypergeometric",
                      use_case = use_case,
                      annotation = annotation,
                      diffRanking = ranking_F - ranking_W,
                      interaction = paste0(method,"-",annotation))
})) 
colnames(fullContent)[1] <- "annotation_id"




# ann_info <- vroom("data/annotation_info_table.tsv") %>%
#   mutate(term = ifelse(nchar(term) > 30, paste0(substr(term,1,30),"..."),term)) %>%
#   as.data.frame()
ann_info <- vroom("data/annotation_info_table.tsv")
ann_info <- rbind(ann_info,c("GO:0042493","response to xenobiotic stimulus"))

annotsStats <- vroom("data/observations/TFsNtargetsPerAnnotationDistribution.tsv") 
colnames(annotsStats)[2] <- "annotation_id"
# dorothea_info <- vroom("data/dorothea.tsv") %>%
#   as.data.frame() %>%
#   filter(org == 9606, confidence <= "C") %>%
#   select(tf,target,confidence) %>%
#   group_by(target) %>%
#   count()

# ann_tables <- rbindlist(lapply(annotations, function(annotation){
#   annotation_table <- gsub(" ","_",annotation)
#   ann_table <- vroom(glue("data/{annotation_table}.tsv")) %>%
#     as.data.frame() %>%
#     filter(organism == 9606, annotation_id %in% ann_info$annotation_id) %>%
#     select(-organism)
#   
#   mean_vals <- rbindlist(mclapply(unique(ann_table$annotation_id), function(ann){
#     filteredTable <- ann_table %>% filter(annotation_id == ann) %>% inner_join(dorothea_info, by = c("symbol" = "target"))
#     data.frame(annotation_id = ann, average = mean(filteredTable$n))
#   }, mc.cores = 10))
#   
#   ann_table <- ann_table %>% inner_join(mean_vals) %>% select(annotation_id,average) %>% unique()
#   ann_table$annotation <- annotation
#   return(ann_table)
# }))

myannotation <- annotations[1]; use_case <- unique(fullContent$use_case)[1]
totalAnnotationsPerDB <- table(fullContent$annotation[fullContent$use_case == use_case]) # == != use_case

topSelected <- 15
subSelection <- fullContent[fullContent$ranking_F %in% 1:topSelected | fullContent$ranking_W %in% 1:topSelected,]
fullInfo <- merge(annotsStats,merge(subSelection,ann_info))

write.table(fullInfo[fullInfo$use_case != use_case & fullInfo$annotation == annotations[2],],"use_case/use_casePlotDataFrame.tsv",sep = "\t",quote = F,row.names = F)

for (use_case in unique(subSelection$use_case)){
  plotList <- lapply(annotations, function(myannotation){
    base_color <- colors_blind[names(colors_blind) %in% myannotation]
    dark_color <- darken(base_color, 0.4)
    light_color <- lighten(base_color,0.4)
    useCaseSpecific <- fullInfo[fullInfo$use_case == use_case & fullInfo$annotation ==  myannotation,]
    useCaseSpecific <- useCaseSpecific[order(useCaseSpecific$diffRanking),]
    useCaseSpecific$term <- factor(useCaseSpecific$term, levels = unique(useCaseSpecific$term))
    useCaseSpecific$annotation_id <- factor(useCaseSpecific$annotation_id, levels = unique(useCaseSpecific$annotation_id))
    #proportionalRank <- diffRanking / totalAnnotationsPerDB[[myannotation]] * 100
    
      ggplot(useCaseSpecific, aes(x = diffRanking, y = annotation_id))+
      geom_col(aes(fill = averageTargetsPerTFInAnnotation),
               width = 0.9, position = "dodge")+
      scale_fill_gradient(high = dark_color, low = light_color)+
      theme(legend.key.size = unit(0.5,"cm"),
            legend.text = element_text(size = 4),
            legend.title = element_text(size = 5),
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(size = 4),
            axis.title.x = element_text(size = 5),
            axis.ticks = element_blank())+
      xlab("Ranking Fisher's Exact Test - Ranking Wallenius's Test")+
      geom_text(aes(label = swr(term,40), hjust = diffRanking > 0), x = 0, size = 1.5, stat = "unique",
                lineheight = 0.7)+
      scale_x_symmetric()+
      labs(fill = "Average TFs\nregulated by\nannotated\ntarget genes")
  })
  
  ggsave(paste0("use_case/",use_case,"_UseCase.tiff"), 
         plot = plot_grid(plotlist = plotList, nrow = 2),
         units = "cm",height = 15, width = 20,
         dpi = 500, compression = "lzw", bg = "white")
  
  ggsave(paste0("use_case/",use_case,"_UseCase.png"), 
         plot = plot_grid(plotlist = plotList, nrow = 2),
         units = "cm",height = 15, width = 20,
         dpi = 500, bg = "white")
}

### OLD IMPLEMENTATION WITH OLDER INPUT FORMAT 
## COMMENTED UNTIL DEFINITIVE IS ACCEPTED
# list_files <- list.files("case_of_use",pattern = "*Cancer*",full.names = T)
# 
# fullContent <- rbindlist(lapply(list_files, function(fileName){
#   fileNameCheck <- gsub("GO_BP","GO BP",fileName)
#   basenameFile <- basename(fileNameCheck)
#   annotation <- gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[2])
#   content <- vroom(fileName) %>% filter(method != "TFs_Hypergeom") %>% mutate(annotation = annotation, interaction = paste0(method,"-",annotation))
# })) %>% mutate(method = factor(method, levels = c("Target_NonCentral","Target_Hypergeom")),
#                annotation = factor(annotation, levels = unique(annotation)))
# 
# 
# annotations = c("KEGG","Reactome","GO BP","WikiPathways")
# 
# colors_blind <- c("KEGG" = "#979A61",
#                   "GO BP" = "#3d91e0",
#                   "WikiPathways" = "#D974A0",
#                   "Reactome" = "#f1b620")
# 
# ann_info <- vroom("data/annotation_info_table.tsv") %>%
#   mutate(term = ifelse(nchar(term) > 30, paste0(substr(term,1,30),"..."),term)) %>%
#   as.data.frame()
# 
# dorothea_info <- vroom("data/dorothea.tsv") %>%
#   as.data.frame() %>%
#   filter(org == 9606, confidence <= "C") %>%
#   select(tf,target,confidence) %>%
#   group_by(target) %>%
#   count()
# 
# ann_tables <- rbindlist(lapply(annotations, function(annotation){
#   annotation_table <- gsub(" ","_",annotation)
#   ann_table <- vroom(glue("data/{annotation_table}.tsv")) %>%
#     as.data.frame() %>%
#     filter(organism == 9606, annotation_id %in% ann_info$annotation_id) %>%
#     select(-organism)
#   
#   mean_vals <- rbindlist(mclapply(unique(ann_table$annotation_id), function(ann){
#     filteredTable <- ann_table %>% filter(annotation_id == ann) %>% inner_join(dorothea_info, by = c("symbol" = "target"))
#     data.frame(annotation_id = ann, average = mean(filteredTable$n))
#   }, mc.cores = 10))
#   
#   ann_table <- ann_table %>% inner_join(mean_vals) %>% select(annotation_id,average) %>% unique()
#   ann_table$annotation <- annotation
#   return(ann_table)
# }))
# 
# plotList <- lapply(annotations, function(annotation){
#   base_color <- colors_blind[names(colors_blind) %in% annotation]
#   dark_color <- darken(base_color, 0.3)
#   light_color <- lighten(base_color,0.3)
#   rankingDifferBar <- fullContent %>% select(annotation_id, method, annotation, ranking) %>% unique() %>% filter(annotation == !!annotation) %>%
#     pivot_wider(names_from = method, values_from = ranking) %>%
#     inner_join(ann_tables) %>%
#     mutate(diffRanking = Target_Hypergeom - Target_NonCentral) %>%
#     filter(Target_NonCentral <= 15 | Target_Hypergeom <= 15) %>%
#     left_join(ann_info) %>%
#     na.omit() %>% arrange(diffRanking) %>% mutate(term = factor(term, levels = unique(term)),
#                                                   annotation_id = factor(annotation_id, levels = unique(annotation_id)),
#                                                   xpos = ifelse(diffRanking < 0, -0.5, 0.5))
#   
#   rankingDifferBar <- rankingDifferBar %>% mutate(prop = diffRanking / length(unique(fullContent %>% filter(annotation == !!annotation) %>% pull(annotation_id))))
#   
#   ggplot(rankingDifferBar, aes(x = diffRanking, y = annotation_id))+
#     geom_col(aes(fill = average))+
#     scale_fill_gradient(high = dark_color, low = light_color)+
#     theme(legend.key.size = unit(0.3,"cm"),
#           legend.text = element_text(size = 4),
#           legend.title = element_text(size = 5),
#           axis.text.y = element_blank(), axis.title.y = element_blank(),
#           panel.background = element_blank(),
#           axis.text.x = element_text(size = 4),
#           axis.title.x = element_text(size = 5),
#           axis.ticks = element_blank())+
#     xlab("Ranking Fisher's Exact Test - Ranking Wallenius's Test")+
#     geom_text(aes(label = term, hjust = diffRanking > 0), x = 0, size = 1.5)+
#     scale_x_symmetric()+
#     labs(fill = "Average TFs\nregulated by\nannotated\ntarget genes")
# })
# 
# fig_6 <- plot_grid(plotlist = plotList, nrow = 2)
#   
# fig5 <- plot_grid(fig_5, fig_6, nrow = 2, labels = c("A","B"),label_fontfamily = "serif", label_size = 10)
# 
# ggsave("Figures/Figure 5.tiff", plot = fig5,
#        units = "cm",height = 25, width = 20,dpi = 500,
#        compression = "lzw", bg = "white")

