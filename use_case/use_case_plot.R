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

currentScriptDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file.path(currentScriptDir,".."))
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")
annotations <- names(colors_blind)


list_files <- list.files("use_case/results",full.names = T)

fullContent <- rbindlist(lapply(list_files, function(fileName){
  fileNameCheck <- gsub("GO_BP","GO BP",fileName)
  basenameFile <- basename(fileNameCheck)
  use_case <- gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[1])
  annotation <- gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[2])
  method <- gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[3])
  content <- vroom(fileName) 
  names(content)[2] <- "pvalue"
  content <- content %>% mutate(method = ifelse(method == "targetsF", "Hypergeometric", "Non-Central Hypergeometric"),
                                use_case = use_case,
                                pvalueAdj = p.adjust(pvalue),
                                ranking = rank(pvalueAdj),
                                annotation = annotation,
                                interaction = paste0(method,"-",annotation))
})) 
fullContent <- fullContent  %>% mutate(annotation = factor(annotation, levels = annotations))

allAnnotationsIDs <- unique(fullContent$term)
for (use_case in unique(fullContent$use_case)){
  for (method in unique(fullContent$method)){
    fullContent[fullContent$use_case == use_case]
  }
}

fullContent_dcasted <- dcast(fullContent, use_case~method)



fullContent_F <- fullContent[fullContent$method == "Hypergeometric", ]
fullContent_W <- fullContent[fullContent$method == "Non-Central Hypergeometric", ]

fullContent_F[match(fullContent_F$term,allAnnotationsIDs)]
fullContent_W[match(fullContent_W$term,allAnnotationsIDs)]



fullContent[fullContent$term %in% unique(fullContent$term)]

table(fullContent$ranking)

ann_info <- vroom("data/annotation_info_table.tsv") %>%
  mutate(term = ifelse(nchar(term) > 30, paste0(substr(term,1,30),"..."),term)) %>%
  as.data.frame()

vroom("data/observations/TFsNtargetsPerAnnotationDistribution.tsv") %>%
  select(annotation, )

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

plotList <- lapply(annotations, function(annotation){
  base_color <- colors_blind[names(colors_blind) %in% annotation]
  dark_color <- darken(base_color, 0.3)
  light_color <- lighten(base_color,0.3)
  rankingDifferBar <- fullContent %>% select(annotation_id, method, annotation, ranking) %>% unique() %>% filter(annotation == !!annotation) %>%
    pivot_wider(names_from = method, values_from = ranking) %>%
    inner_join(ann_tables) %>%
    mutate(diffRanking = Target_Hypergeom - Target_NonCentral) %>%
    filter(Target_NonCentral <= 15 | Target_Hypergeom <= 15) %>%
    left_join(ann_info) %>%
    na.omit() %>% arrange(diffRanking) %>% mutate(term = factor(term, levels = unique(term)),
                                                  annotation_id = factor(annotation_id, levels = unique(annotation_id)),
                                                  xpos = ifelse(diffRanking < 0, -0.5, 0.5))
  
  rankingDifferBar <- rankingDifferBar %>% mutate(prop = diffRanking / length(unique(fullContent %>% filter(annotation == !!annotation) %>% pull(annotation_id))))
  
  ggplot(rankingDifferBar, aes(x = diffRanking, y = annotation_id))+
    geom_col(aes(fill = average))+
    scale_fill_gradient(high = dark_color, low = light_color)+
    theme(legend.key.size = unit(0.3,"cm"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 4),
          axis.title.x = element_text(size = 5),
          axis.ticks = element_blank())+
    xlab("Ranking Fisher's Exact Test - Ranking Wallenius's Test")+
    geom_text(aes(label = term, hjust = diffRanking > 0), x = 0, size = 1.5)+
    scale_x_symmetric()+
    labs(fill = "Average TFs\nregulated by\nannotated\ntarget genes")
})

fig_5 <- plot_grid(plotlist = plotList, nrow = 2)

list_files <- list.files("case_of_use",pattern = "*Cancer*",full.names = T)

fullContent <- rbindlist(lapply(list_files, function(fileName){
  fileNameCheck <- gsub("GO_BP","GO BP",fileName)
  basenameFile <- basename(fileNameCheck)
  annotation <- gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[2])
  content <- vroom(fileName) %>% filter(method != "TFs_Hypergeom") %>% mutate(annotation = annotation, interaction = paste0(method,"-",annotation))
})) %>% mutate(method = factor(method, levels = c("Target_NonCentral","Target_Hypergeom")),
               annotation = factor(annotation, levels = unique(annotation)))


annotations = c("KEGG","Reactome","GO BP","WikiPathways")

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

ann_info <- vroom("data/annotation_info_table.tsv") %>%
  mutate(term = ifelse(nchar(term) > 30, paste0(substr(term,1,30),"..."),term)) %>%
  as.data.frame()

dorothea_info <- vroom("data/dorothea.tsv") %>%
  as.data.frame() %>%
  filter(org == 9606, confidence <= "C") %>%
  select(tf,target,confidence) %>%
  group_by(target) %>%
  count()

ann_tables <- rbindlist(lapply(annotations, function(annotation){
  annotation_table <- gsub(" ","_",annotation)
  ann_table <- vroom(glue("data/{annotation_table}.tsv")) %>%
    as.data.frame() %>%
    filter(organism == 9606, annotation_id %in% ann_info$annotation_id) %>%
    select(-organism)
  
  mean_vals <- rbindlist(mclapply(unique(ann_table$annotation_id), function(ann){
    filteredTable <- ann_table %>% filter(annotation_id == ann) %>% inner_join(dorothea_info, by = c("symbol" = "target"))
    data.frame(annotation_id = ann, average = mean(filteredTable$n))
  }, mc.cores = 10))
  
  ann_table <- ann_table %>% inner_join(mean_vals) %>% select(annotation_id,average) %>% unique()
  ann_table$annotation <- annotation
  return(ann_table)
}))

plotList <- lapply(annotations, function(annotation){
  base_color <- colors_blind[names(colors_blind) %in% annotation]
  dark_color <- darken(base_color, 0.3)
  light_color <- lighten(base_color,0.3)
  rankingDifferBar <- fullContent %>% select(annotation_id, method, annotation, ranking) %>% unique() %>% filter(annotation == !!annotation) %>%
    pivot_wider(names_from = method, values_from = ranking) %>%
    inner_join(ann_tables) %>%
    mutate(diffRanking = Target_Hypergeom - Target_NonCentral) %>%
    filter(Target_NonCentral <= 15 | Target_Hypergeom <= 15) %>%
    left_join(ann_info) %>%
    na.omit() %>% arrange(diffRanking) %>% mutate(term = factor(term, levels = unique(term)),
                                                  annotation_id = factor(annotation_id, levels = unique(annotation_id)),
                                                  xpos = ifelse(diffRanking < 0, -0.5, 0.5))
  
  rankingDifferBar <- rankingDifferBar %>% mutate(prop = diffRanking / length(unique(fullContent %>% filter(annotation == !!annotation) %>% pull(annotation_id))))
  
  ggplot(rankingDifferBar, aes(x = diffRanking, y = annotation_id))+
    geom_col(aes(fill = average))+
    scale_fill_gradient(high = dark_color, low = light_color)+
    theme(legend.key.size = unit(0.3,"cm"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 4),
          axis.title.x = element_text(size = 5),
          axis.ticks = element_blank())+
    xlab("Ranking Fisher's Exact Test - Ranking Wallenius's Test")+
    geom_text(aes(label = term, hjust = diffRanking > 0), x = 0, size = 1.5)+
    scale_x_symmetric()+
    labs(fill = "Average TFs\nregulated by\nannotated\ntarget genes")
})

fig_6 <- plot_grid(plotlist = plotList, nrow = 2)
  
fig5 <- plot_grid(fig_5, fig_6, nrow = 2, labels = c("A","B"),label_fontfamily = "serif", label_size = 10)

ggsave("Figures/Figure 5.tiff", plot = fig5,
       units = "cm",height = 25, width = 20,dpi = 500,
       compression = "lzw", bg = "white")

