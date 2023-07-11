library(vroom)
library(glue)
library(tidyverse)
library(parallel)
library(xlsx)
library(ggplot2)
library(sdamr)
library(ggrepel)
library(cowplot)
library(stringr)
library(data.table)
library(colorspace)
library(ggpubr)
library(lemon)
library(gtable)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

list_files <- list.files("case_of_use",pattern = "*SLE*",full.names = T)

fullContent <- rbindlist(lapply(list_files, function(fileName){
  fileNameCheck <- gsub("GO_BP","GO BP",fileName)
  basenameFile <- basename(fileNameCheck)
  annotation <- gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[2])
  content <- vroom(fileName) %>% filter(typeOfAnalysis != "TFs_Hypergeom") %>% mutate(annotation = annotation, interaction = paste0(typeOfAnalysis,"-",annotation))
})) %>% mutate(typeOfAnalysis = factor(typeOfAnalysis, levels = c("Target_NonCentral","Target_Hypergeom")),
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
  rankingDifferBar <- fullContent %>% select(annotation_id, typeOfAnalysis, annotation, ranking) %>% unique() %>% filter(annotation == !!annotation) %>%
    pivot_wider(names_from = typeOfAnalysis, values_from = ranking) %>%
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
  content <- vroom(fileName) %>% filter(typeOfAnalysis != "TFs_Hypergeom") %>% mutate(annotation = annotation, interaction = paste0(typeOfAnalysis,"-",annotation))
})) %>% mutate(typeOfAnalysis = factor(typeOfAnalysis, levels = c("Target_NonCentral","Target_Hypergeom")),
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
  rankingDifferBar <- fullContent %>% select(annotation_id, typeOfAnalysis, annotation, ranking) %>% unique() %>% filter(annotation == !!annotation) %>%
    pivot_wider(names_from = typeOfAnalysis, values_from = ranking) %>%
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



# fig_6 <- ggplot(rankingDiffer, aes(y = diffRanking, x = average))+
#   facet_wrap(~annotation, scales = "free")+
#   geom_point(size = 0.5, aes(color = annotation))+
#   scale_color_manual(values = colors_blind)+
#   # geom_smooth(linewidth = 0.3, aes(color = annotation))+
#   theme(panel.background = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.background = element_rect(fill = "white"),
#         axis.text.x = element_text(size = 5),
#         title = element_text(face = "bold"),
#         strip.text = element_text(size = 5,face = "bold"),
#         axis.title = element_text(size = 6),
#         axis.text.y = element_text(size = 5),
#         axis.line.x = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
#         legend.margin=margin(0,8,0,0),
#         legend.box.margin=margin(-10,0,-10,-10),
#         legend.key = element_rect(fill = "white"),
#         legend.text = element_text(size = 5),
#         legend.title = element_text(size = 6),
#         legend.position = "none")+
#   xlab("Average TFs regulated by annotated target genes")+
#   ylab("Ranking Fisher's Exact Test - Ranking Wallenius's Test")+
#   geom_text_repel(aes(label = term), size = 1.5,segment.size = 0.1)
# 
# 
# ggplot(rankingDiffer, aes(x = diffRanking, y = term))+
#   geom_col(aes(fill = average))+
#   facet_wrap(~annotation, scales = "free")+
#   theme(axis.text.y = element_blank(),
#         axis.title = element_blank())+
#   geom_text(aes(label = term))+
#   scale_fill_gradient(low = "pink", high = "darkred")
fig_6 <- plot_grid(plotlist = plotList, nrow = 2)
  

fig5 <- plot_grid(fig_5, fig_6, nrow = 2, labels = c("A","B"),label_fontfamily = "serif", label_size = 10)

ggsave("Figures/Figure 5.tiff", plot = fig5,
       units = "cm",height = 25, width = 20,dpi = 500,
       compression = "lzw", bg = "white")

