library(vroom)
library(glue)
library(dplyr)
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
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

### Prepare data

annotations = c("KEGG","Reactome","GO BP","WikiPathways")

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

ann_info <- vroom("data/annotation_info_table.tsv") %>%
  as.data.frame()

ann_tables <- do.call("rbind",mclapply(annotations, function(annotation){
  annotation_table <- gsub(" ","_",annotation)
  ann_table <- vroom(glue("data/{annotation_table}.tsv")) %>%
    as.data.frame() %>%
    filter(organism == 9606, annotation_id %in% ann_info$annotation_id) %>%
    select(-organism)
  ann_table$annotation <- annotation
  return(ann_table)
}))

dorothea_info <- vroom("data/dorothea.tsv") %>%
  as.data.frame() %>%
  filter(org == 9606) %>%
  select(tf,target,confidence)

trrust <- vroom("trrust_rawdata.human.tsv",col_names = F)[,1:2] %>%
  as.data.frame()
colnames(trrust) <- c("tf","target")
trrust$confidence <- "TRRUST v2"

regulons <- rbind(dorothea_info,trrust)

################################ Plot Gene-Term Pairs (%). Figure 2.tiff

regulons_coverage <- do.call("rbind", lapply(annotations, function(ann){
  ann_table <- ann_tables[ann_tables$annotation == ann,]
  
  ## Number Gene-Term Pairs
  n <- nrow(unique(ann_table))
  
  confidence_level <- do.call("rbind",mclapply(c("TRRUST v2","A","B","C","D","E"), function(conf){
    if (conf != "TRRUST v2"){
      regulons_filt <- regulons[regulons$confidence <= conf,] 
    } else{
      regulons_filt <- regulons[regulons$confidence == conf,] 
    }
    target_genes <- unique(regulons_filt$target)
    target_genes <- ann_table[ann_table$symbol %in% target_genes,]
    n_genes <- length(unique(target_genes$symbol))
    target_genes <- round((nrow(unique(target_genes))/n)*100,2)
    if (confidence != "TRRUST v2"){
      confidence = paste0("DoRothEA ",confidence)
    }
    target_genes <- data.frame("coverage" = target_genes,
                               "n_genes" = n_genes,
                               "confidence" = confidence,
                               "type" = "Target Genes")
  }))
  
  tfs_dorothea <- regulons %>%
    filter(confidence <= "E") %>%
    select(tf) %>%
    unique() %>%
    pull()
  
  tfs_dorothea <- ann_table[ann_table$symbol %in% tfs_dorothea,]
  n_genes <- length(unique(tfs_dorothea$symbol))
  tfs_dorothea <- round((nrow(unique(tfs_dorothea)) / n) * 100,2)
  tfs_dorothea <- data.frame("coverage" = tfs_dorothea,
                             "n_genes" = n_genes,
                             "confidence" = "DoRothEA TFs",
                             "type" = "TFs")
  
  tfs_trrust <- regulons %>%
    filter(confidence == "TRRUST v2") %>%
    select(tf) %>%
    unique() %>%
    pull()
  
  tfs_trrust <- ann_table[ann_table$symbol %in% tfs_trrust,]
  n_genes <- length(unique(tfs_trrust$symbol))
  tfs_trrust <- round((nrow(unique(tfs_trrust)) / n) * 100,2)
  tfs_trrust <- data.frame("coverage" = tfs_trrust,
                             "n_genes" = n_genes,
                             "confidence" = "TRRUST v2 TFs",
                             "type" = "TFs")
  
  coverage <- rbind(tfs_trrust, tfs_dorothea, confidence_level)
  coverage$annotation = annotation
  
  to_file <- coverage[,c("confidence","coverage","n_genes","type")]
  colnames(to_file) <- c("TFs / Target confidence", "% coverage annotation", "Number of genes","Type")
  
  if (annotation == "KEGG"){
    invisible(file.remove("Tables/Supplementary Table 1.xlsx"))
    write.xlsx(to_file, file = "Tables/Supplementary Table 1.xlsx",
               append = FALSE,
               sheetName = annotation,
               row.names = F)
  } else{
    write.xlsx(to_file, file = "Tables/Supplementary Table 1.xlsx",
               append = T,
               sheetName = annotation,
               row.names = F)
  }
  return(coverage)
}))

regulons_coverage$annotation <- factor(regulons_coverage$annotation, levels = annotations)
regulons_coverage$type <- factor(regulons_coverage$type, levels = c("TFs","Target Genes"))
regulons_coverage$confidence <- factor(regulons_coverage$confidence, levels = unique(regulons_coverage$confidence))

gg_2 <- ggplot(regulons_coverage, aes(y = coverage, x = confidence, colour = annotation))+
  facet_wrap(~type, scales = "free_x")+
  geom_line(aes(group = annotation), size = 0.2)+
  geom_point(size = 2.5,colour = "black",alpha = 0.5)+
  geom_point(size = 2) +
  scale_color_manual(values = colors_blind)+
  theme(panel.background = element_blank(),
        strip.text = element_text(size = 5),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 5),
        axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 4),
        legend.key = element_rect(fill = "white"),
        legend.key.size  = unit(0.1, 'cm'),
        # legend.title = element_text(size = 6),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.position = "bottom")+
  guides(colour = guide_legend(override.aes = list(size = 1),
                               title = "Annotation"))+
  ylab("Annotation coverage (%)")+
  xlab("")

ggsave("Figures/Figure 2.tiff", plot = gg_2,units = "cm",height = 8, width = 8,dpi = 500,
       compression = "lzw",bg = "white")

########################################################

### Use TFs directly Supplementary Figure 1 & Supplementary Table 2

content <- lapply(annotations, function(annotation){
  message(glue("Reading files from {annotation}"))
  nsize = c(2,3,4,5,10,20,30)

  content <- mclapply(nsize, function(n){
    result_file = paste0("null_simulation_tfs_rank/",gsub(" ","_",annotation),"_",n,".tsv")
    full_content <- vroom(result_file)

    content = full_content %>%
      group_by(term) %>%
      summarise(median_rank_hypergeom = median(rank_hyper),
                proportion_hypergeom = sum(pval_hyper < 0.05) / 1000)

    content$n = n
    full_content$n = n
    content <- list(content,full_content)
    return(content)
  })

  message(glue("Joining files from {annotation}"))
  full_content <- rbindlist(sapply(content,`[`,2))
  content <-  rbindlist(sapply(content,`[`,1))

  content$n <- factor(content$n, levels = unique(content$n))
  content$annotation = annotation
  full_content$n <- factor(full_content$n, levels = unique(full_content$n))
  full_content$annotation = annotation
  content <- left_join(content, ann_info, by = c("term" = "annotation_id"))
  full_content <- left_join(full_content, ann_info, by = c("term" = "annotation_id"))

  to_file <- content[order(content$proportion_hypergeom,decreasing = T),c("term","term.y","n","proportion_hypergeom","median_rank_hypergeom")]
  colnames(to_file) <- c("Annotation Id","Term","Number of TFs","Proportion of times significant (p < 0.05)","Median Ranking")

  options(java.parameters = "-Xmx8000m")
  invisible(gc())
  message(glue("Writing table for {annotation}"))

  write.xlsx(to_file,file = glue("Tables/Supplementary Table 2 {annotation}.xlsx"),
             append = FALSE,
             sheetName = annotation,
             row.names = F)

  content = list(content, full_content)
  message(glue("Closing {annotation}"))
  return(content)
})

full_content <- rbindlist(sapply(content,`[`,2))
content <-  rbindlist(sapply(content,`[`,1))

content$annotation <- factor(content$annotation, levels = annotations)

ggsup_1_v1 <- ggplot(content, aes(x = n, y = median_rank_hypergeom, color = annotation))+
  # geom_boxjitter(position = position_dodge(width = 1), jitter.size = 0.05, outlier.shape = NA, lwd = 0.2)+
  geom_point(position = position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7),
             size = 0.2)+
  geom_boxplot(outlier.shape = NA, lwd = 0.4, alpha = 0.5)+
  facet_wrap(~annotation, scales = "free", nrow = 1)+
  scale_color_manual(values = colors_blind)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_blank())+
  xlab("Number of TFs")+
  ylab("Median Ranking")

ggsup_1_v2 <- ggplot(content, aes(x = n, y = proportion_hypergeom, color = annotation))+
  # geom_boxjitter(position = position_dodge(width = 1), jitter.size = 0.05, outlier.shape = NA, lwd = 0.2)+
  geom_point(position = position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7),
             size = 0.2)+
  geom_boxplot(outlier.shape = NA, lwd = 0.4, alpha = 0.5)+
  facet_wrap(~annotation, scales = "free", nrow = 1)+
  scale_color_manual(values = colors_blind)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_blank())+
  xlab("Number of TFs")+
  ylab("Proportion of times significant (p < 0.05)")


content_top = content %>%
  filter(n == 30) %>%
  arrange(median_rank_hypergeom) %>%
  group_by(annotation) %>%
  do(head(.,3))

content_top <- full_content[full_content$term %in% content_top$term,]
content_top$annotation <- factor(content_top$annotation, levels = annotations)
content_top <- content_top[order(content_top$annotation),]
content_top$term.y <- ifelse(content_top$term != "WP706", content_top$term.y,"Sudden infant death syndrome (SIDS)\n susceptibility pathways")
content_top$term.y <- paste0(content_top$term.y," (",content_top$term,")")
content_top$term.y <- ifelse(content_top$term %in% c("R-HSA-6785807","WP3611","WP3617") | content_top$annotation == "GO BP", gsub(" \\(","\n(",content_top$term.y),content_top$term.y)
content_top$term.y <- factor(content_top$term.y, levels = unique(content_top$term.y))

ggsup_2_v1 <- ggplot(content_top, aes(y = rank_hyper, x = n, color = annotation))+
  facet_wrap(~term.y, nrow = 4, scales = "free")+
  geom_point(position = position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7),
             size = 0.3)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.4)+
  scale_color_manual(values = colors_blind)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.position = "none")+
  xlab("Number of TFs")+
  ylab("Rank")


content_top = content %>%
  filter(n == 30) %>%
  arrange(desc(proportion_hypergeom)) %>%
  group_by(annotation) %>%
  do(head(.,3))

content_top <- full_content[full_content$term %in% content_top$term,]
content_top$annotation <- factor(content_top$annotation, levels = annotations)
content_top <- content_top[order(content_top$annotation),]
content_top$term.y <- ifelse(content_top$term != "WP706", content_top$term.y,"Sudden infant death syndrome (SIDS)\n susceptibility pathways")
content_top$term.y <- paste0(content_top$term.y," (",content_top$term,")")
content_top$term.y <- ifelse(content_top$term %in% c("R-HSA-6785807","WP3611","WP3617") | content_top$annotation == "GO BP", gsub(" \\(","\n(",content_top$term.y),content_top$term.y)
content_top$term.y <- factor(content_top$term.y, levels = unique(content_top$term.y))


ggsup_2_v2 <- ggplot(content_top, aes(y = -log(pval_hyper), x = n, color = annotation))+
  facet_wrap(~term.y, nrow = 4, scales = "free")+
  geom_point(position = position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7),
             size = 0.3)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.4)+
  scale_color_manual(values = colors_blind)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.position = "none")+
  xlab("Number of TFs")+
  ylab("- log (pvalue)")

gg_legend <- get_legend(ggsup_1_v2)
ggsup_1_v2 <- ggsup_1_v2 + theme(legend.position = "none")
ggsup_v2 <- plot_grid(ggsup_1_v2, ggsup_2_v2, nrow = 2, rel_heights = c(0.8,1.2),labels = c('A', 'B'),
                      label_fontfamily = "serif", label_size = 10,
                      label_fontface = "plain")
ggsup_v2 <- plot_grid(ggsup_v2, gg_legend, rel_widths = c(0.9,0.2))
ggsave("Figures/Supplementary Figure 1.tiff", plot = ggsup_v2,
       units = "cm",height = 18, width = 20,dpi = 500,
       compression = "lzw", bg = "white")


gg_legend <- get_legend(ggsup_1_v1)
ggsup_1_v1 <- ggsup_1_v1 + theme(legend.position = "none")
ggsup_v1 <- plot_grid(ggsup_1_v1, ggsup_2_v1, nrow = 2, rel_heights = c(0.8,1.2),labels = c('A', 'B'),
                      label_fontfamily = "serif", label_size = 10,
                      label_fontface = "plain")
ggsup_v1 <- plot_grid(ggsup_v1, gg_legend, rel_widths = c(0.9,0.2))
ggsave("Figures/Supplementary Figure 2.tiff", plot = ggsup_v1,
       units = "cm",height = 18, width = 20,dpi = 500,
       compression = "lzw", bg = "white")

#######################################################


dorothea_info <- dorothea_info %>%
  filter(confidence <= "C")

count_genes <- do.call("rbind",lapply(annotations, function(annotation){
  data_content <- ann_tables[ann_tables$annotation == annotation,]
  target_genes <- unique(dorothea_info$target)
  target_genes <- data_content[data_content$symbol %in% target_genes,]
  target_genes <- unique(target_genes$symbol)
  dorothea_filt <- dorothea_info[dorothea_info$target %in% target_genes,] %>%
    group_by(target) %>%
    summarise(tfs = n()) %>%
    mutate(annotation = annotation)

  to_file <- dorothea_filt[order(dorothea_filt$tfs, decreasing = T),c("target","tfs")] %>%
    as.data.frame()
  colnames(to_file) <- c("Target Gene","Number of TFs")

  options(java.parameters = "-Xmx8000m")
  message(glue("Writing table for {annotation}"))
  if (annotation == "KEGG"){
    invisible(file.remove("Tables/Supplementary Table 3.xlsx"))
    write.xlsx(to_file, file = "Tables/Supplementary Table 3.xlsx",
               append = FALSE,
               sheetName = annotation,
               row.names = F)
  } else{
    write.xlsx(to_file, file = "Tables/Supplementary Table 3.xlsx",
               append = T,
               sheetName = annotation,
               row.names = F)
  }
  return(dorothea_filt)
}))

count_genes_top <- count_genes %>%
  arrange(desc(tfs)) %>%
  group_by(annotation) %>%
  do(head(.,20))

count_genes_top_list <- unique(count_genes_top$target)
count_genes$annotation <- factor(count_genes$annotation, levels = rev(unique(count_genes$annotation)))

gg3_1 <- ggplot(count_genes, aes(annotation,y = tfs,fill=annotation))+
  geom_flat_violin(position = position_nudge(x = .05, y = 0),
                   alpha = .8, adjust = 2) +
  geom_point(aes(y = tfs, color = annotation), position = position_jitternudge(nudge.x = -0.3,jitter.width = 0.2,seed = 123), size = .5, alpha = 0.8) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = -.3, y = 0))+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "none",
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 6),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))+
  scale_color_manual(values = colors_blind) +
  scale_fill_manual(values = colors_blind) +
  ylab("TFs that regulate each gene")+
  geom_text_repel(data = count_genes[count_genes$target %in% count_genes_top_list,],aes(annotation,tfs,label=target),
                  position = position_jitternudge(nudge.x = -0.3,jitter.width = 0.1,seed = 123),size = 1.5,
                  min.segment.length = 0.3, max.overlaps = 20)+
  coord_flip()

ann_table <- do.call("rbind",lapply(annotations, function(annotation){
  count_genes_top <- count_genes_top[count_genes_top$annotation == annotation,]
  ann_table <- ann_tables[ann_tables$annotation == annotation & ann_tables$symbol %in% count_genes_top$target,]
  ann_table_show <- as.data.frame(table(ann_table$annotation_id)) %>%
    arrange(desc(Freq)) %>%
    do(head(.,10))
  ann_table_show$annotation_id <- as.character(ann_table_show$Var1)
  ann_table_show$annotation <- annotation
  ann_table_show <- merge(ann_table_show,ann_info,by = "annotation_id") %>%
    select(annotation_id,annotation,term,Freq)

  to_file <- ann_table_show[order(ann_table_show$Freq,decreasing = T),c("annotation_id","term","Freq")] %>%
    as.data.frame()
  colnames(to_file) <- c("Annotation Id","Term","Number of top genes implicated in gene set")

  options(java.parameters = "-Xmx8000m")
  message(glue("Writing table for {annotation}"))
  if (annotation == "KEGG"){
    invisible(file.remove("Tables/Supplementary Table 4.xlsx"))
    write.xlsx(to_file, file = "Tables/Supplementary Table 4.xlsx",
               append = FALSE,
               sheetName = annotation,
               row.names = F)
  } else{
    write.xlsx(to_file, file = "Tables/Supplementary Table 4.xlsx",
               append = T,
               sheetName = annotation,
               row.names = F)
  }
  return(ann_table_show)
}))

ann_table <- ann_table[order(ann_table$Freq,decreasing = T),]
ann_table$annotation_id <- factor(ann_table$annotation_id,levels = unique(ann_table$annotation_id))

gg3_2 <- ggplot(ann_table, aes(x = Freq, y = annotation_id, fill = annotation, label = term)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colors_blind)+
  scale_x_continuous(expand = c(0,0),limits = c(0,max(ann_table$Freq)+2))+
  geom_text(size = 2,x = 1,hjust = 0)+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key = element_rect(fill = "white"),
        legend.position = "bottom",
        legend.key.size = unit(0.3,"cm"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        plot.margin = margin(0.2, 0.2, 0.2, 0.5, "cm"))+
  labs(fill = "Annotation")+
  xlab("Target genes related to each term")

gg_legend <- get_legend(gg3_2)

gg3_2 <- gg3_2 + theme(legend.position = "none")

gg <- plot_grid(gg3_1,gg3_2, labels = c("A","B"),
                label_fontfamily = "serif", label_size = 10,
                label_fontface = "plain")

gg <- plot_grid(gg,gg_legend, rel_heights = c(0.9,0.1), nrow = 2)

ggsave("Figures/Figure 3.tiff", plot = gg,units = "cm",height = 12, width = 20,dpi = 500,
       compression = "lzw",bg = "white")
# 
# 
# ##########################################################################
# 
content <- lapply(annotations, function(annotation){
  message(glue("Reading files from {annotation}"))
  nsize = c(2,3,4,5,10,20,30)

  content <- mclapply(nsize, function(n){
    result_file = paste0("null_simulation_wallenius_rank/",gsub(" ","_",annotation),"_",n,".tsv")
    full_content <- vroom(result_file)

    content = full_content %>%
      group_by(term) %>%
      summarise(median_rank_hypergeom = median(rank_hyper),
                median_rank_noncentral = median(rank_noncentral),
                proportion_hypergeom = sum(hyper < 0.05) / 1000,
                proportion_noncentral = sum(noncentral < 0.05) / 1000)

    content <- rbindlist(lapply(c("hypergeom","noncentral"), function(method){
      columns <- grepl(method,colnames(content))
      columns <- colnames(content)[columns]
      content <- content[,c("term",columns)]
      median_rank <- content[,str_detect(colnames(content),"median")]
      colnames(median_rank) <- "median_rank"
      median_rank <- median_rank$median_rank
      proportion <- content[,str_detect(colnames(content),"proportion")]
      colnames(proportion) <- "proportion"
      proportion <- proportion$proportion
      terms <- content$term
      content <- data.frame(term = terms, Method = method, median_rank = median_rank, proportion = proportion)
    }))

    content$n = n
    content$Method <- ifelse(content$Method == "hypergeom",
                             "Fisher's Exact Test","Wallenius' Test")

    full_content <- rbindlist(lapply(c("hyper","noncentral"), function(method){
      columns <- grepl(method,colnames(full_content))
      columns <- colnames(full_content)[columns]
      content <- full_content[,c("term",columns)]
      rank_vals <- content[,str_detect(colnames(content),"rank")]
      colnames(rank_vals) <- "rank"
      rank_vals <- rank_vals$rank
      pval <- content[,method]
      colnames(pval) <- "pval"
      pval <- pval$pval
      terms <- content$term
      content <- data.frame(term = terms, Method = method, rank = rank_vals, pval = pval)
    }))

    full_content$n = n
    full_content$Method <- ifelse(full_content$Method == "hyper",
                                  "Fisher's Exact Test","Wallenius' Test")
    content <- list(content,full_content)
    return(content)
  })
  message(glue("Joining files from {annotation}"))
  full_content <- rbindlist(sapply(content,`[`,2))
  content <- rbindlist(sapply(content,`[`,1))

  content$n <- factor(content$n, levels = unique(content$n))
  content$annotation = annotation
  full_content$n <- factor(full_content$n, levels = unique(full_content$n))
  full_content$annotation = annotation
  content <- left_join(content, ann_info, by = c("term" = "annotation_id"))
  full_content <- left_join(full_content, ann_info, by = c("term" = "annotation_id"))

  to_file_hyper <- content[content$Method == "Fisher's Exact Test",c("term","term.y","n","proportion","median_rank")] %>%
    as.data.frame()
  to_file_hyper <- to_file_hyper[order(to_file_hyper$proportion,decreasing = T),]
  colnames(to_file_hyper) <- c("Annotation Id","Term","Number of TFs","Proportion of times significant (p < 0.05)","Median Rank")

  to_file_noncentral <- content[content$Method != "Fisher's Exact Test",c("term","term.y","n","proportion","median_rank")] %>%
    as.data.frame()
  to_file_noncentral <- to_file_noncentral[order(to_file_noncentral$proportion,decreasing = T),]
  colnames(to_file_noncentral) <- c("Annotation Id","Term","Number of TFs","Proportion of times significant (p < 0.05)","Median Rank")

  options(java.parameters = "-Xmx8000m")
  message(glue("Writing table for {annotation}"))
  invisible(gc())
  write.table(to_file_hyper, file = glue("Tables/Supplementary Table 5 Fisher's Exact Test {annotation}.tsv"),
               row.names = F, quote = F, sep = "\t")
  write.table(to_file_noncentral, file = glue("Tables/Supplementary Table 5 Wallenius' Test {annotation}.tsv"),
              row.names = F, quote = F, sep = "\t")

  content = list(content, full_content)
  message(glue("Closing {annotation}"))
  return(content)
})

color_vals <- c()
for (annotation in annotations){
  base_color <- colors_blind[names(colors_blind) %in% annotation]
  dark_color <- darken(base_color, 0.3)
  light_color <- lighten(base_color,0.3)
  color_vals_insert <- c(dark_color, light_color)
  names(color_vals_insert) <- c(paste0(annotation," Fisher's Exact Test"),
                                paste0(annotation," Wallenius' Test"))
  color_vals <- c(color_vals, color_vals_insert)
}

full_content <- rbindlist(sapply(content,`[`,2))
content <- rbindlist(sapply(content,`[`,1))

content$annotation <- factor(content$annotation, levels = annotations)
content$intersection <- paste0(content$annotation, " ", content$Method)
content$intersection <- factor(content$intersection, levels = unique(content$intersection))

full_content$annotation <- factor(full_content$annotation, levels = annotations)
full_content$intersection <- paste0(full_content$annotation, " ", full_content$Method)
full_content$intersection <- factor(full_content$intersection, levels = unique(full_content$intersection))

gg4_1_v1 <- ggplot(content, aes(x = n, y = median_rank, color = intersection))+
  # geom_boxjitter(position = position_dodge(width = 1), jitter.size = 0.05, outlier.shape = NA, lwd = 0.2)+
  geom_point(position = position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7),
             size = 0.2)+
  geom_boxplot(outlier.shape = NA, lwd = 0.4, alpha = 0.5)+
  facet_wrap(~annotation, scales = "free", nrow = 1)+
  scale_color_manual(values = color_vals)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_blank())+
  xlab("Number of TFs")+
  ylab("Median Ranking")

gg4_1_v2 <- ggplot(content, aes(x = n, y = proportion, color = intersection))+
  # geom_boxjitter(position = position_dodge(width = 1), jitter.size = 0.05, outlier.shape = NA, lwd = 0.2)+
  geom_point(position = position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7),
             size = 0.2)+
  geom_boxplot(outlier.shape = NA, lwd = 0.4, alpha = 0.5)+
  facet_wrap(~annotation, scales = "free", nrow = 1)+
  scale_color_manual(values = color_vals)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_blank())+
  xlab("Number of TFs")+
  ylab("Proportion of times significant (p < 0.05)")

content_top = content %>%
  filter(n == 30) %>%
  arrange(median_rank) %>%
  group_by(annotation) %>%
  do(head(.,3))

content_top <- full_content[full_content$term %in% content_top$term,]
content_top$annotation <- factor(content_top$annotation, levels = annotations)
content_top <- content_top[order(content_top$annotation),]
content_top$term.y <- paste0(content_top$term.y," (",content_top$term,")")
content_top$term.y <- ifelse(content_top$term %in% c("R-HSA-6785807","WP3611","WP3617") | content_top$annotation == "GO BP", gsub(" \\(","\n(",content_top$term.y),content_top$term.y)
content_top$term.y <- factor(content_top$term.y, levels = unique(content_top$term.y))

gg4_2_v1 <- ggplot(content_top, aes(y = rank,x = n, color = intersection))+
  facet_wrap(~term.y, nrow = 4, scales = "free")+
  geom_point(position = position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7),
             size = 0.3)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.4)+
  scale_color_manual(values = color_vals)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.position = "none")+
  xlab("Number of TFs")+
  ylab("Rank")


content_top = content %>%
  filter(n == 30) %>%
  arrange(desc(proportion)) %>%
  group_by(annotation) %>%
  do(head(.,3))

content_top <- full_content[full_content$term %in% content_top$term,]
content_top$annotation <- factor(content_top$annotation, levels = annotations)
content_top <- content_top[order(content_top$annotation),]
content_top$term.y <- paste0(content_top$term.y," (",content_top$term,")")
content_top$term.y <- ifelse(content_top$term %in% c("R-HSA-6785807","WP3611") | content_top$annotation == "GO BP", gsub(" \\(","\n(",content_top$term.y),content_top$term.y)
content_top$term.y <- factor(content_top$term.y, levels = unique(content_top$term.y))


gg4_2_v2 <- ggplot(content_top, aes(y = -log(pval), x = n, color = intersection))+
  facet_wrap(~term.y, nrow = 4, scales = "free")+
  geom_point(position = position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7),
             size = 0.3)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.4)+
  scale_color_manual(values = color_vals)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.position = "none")+
  xlab("Number of TFs")+
  ylab("- log (pvalue)")

gg_legend <- get_legend(gg4_1_v2)
gg4_1_v2 <- gg4_1_v2 + theme(legend.position = "none")
gg4_v2 <- plot_grid(gg4_1_v2, gg4_2_v2, nrow = 2, rel_heights = c(0.8,1.2),labels = c('A', 'B'),
                      label_fontfamily = "serif", label_size = 10,
                      label_fontface = "plain")
gg4_v2 <- plot_grid(gg4_v2, gg_legend, rel_widths = c(0.9,0.2))
ggsave("Figures/Figure 4.tiff", plot = gg4_v2,
       units = "cm",height = 18, width = 20,dpi = 500,
       compression = "lzw", bg = "white")


gg_legend <- get_legend(gg4_1_v1)
gg4_1_v1 <- gg4_1_v1 + theme(legend.position = "none")
gg4_v1 <- plot_grid(gg4_1_v1, gg4_2_v1, nrow = 2, rel_heights = c(0.8,1.2),labels = c('A', 'B'),
                      label_fontfamily = "serif", label_size = 10,
                      label_fontface = "plain")
gg4_v1 <- plot_grid(gg4_v1, gg_legend, rel_widths = c(0.9,0.2))
ggsave("Figures/Supplementary Figure 3.tiff", plot = gg4_v1,
       units = "cm",height = 18, width = 20,dpi = 500,
       compression = "lzw", bg = "white")


##### Case of Use

sle_results <- rbindlist(lapply(annotations, function(annotation){
  result_file <- paste0("case_of_use/SLE_",gsub(" ","_",annotation),".tsv")
  result_file <- vroom(result_file) %>%
    as.data.frame()
  
  result_file <- rbindlist(lapply(c("hyper","noncentral"), function(method){
    columns <- grepl(method,colnames(result_file))
    columns <- colnames(result_file)[columns]
    content <- result_file[,c("term",columns)]
    rank_vals <- content[,str_detect(colnames(content),"rank")]
    pval <- content[,method]
    terms <- content$term
    content <- data.frame(term = terms, Method = method, rank = rank_vals, pval = pval)
  }))
  result_file$annotation <- annotation
  result_file$Method <- ifelse(result_file$Method == "hyper",
                               "Fisher's Exact Test","Wallenius' Test")
  result_file <- left_join(result_file, ann_info, by = c("term" = "annotation_id"))
  to_file_hyper <- result_file[result_file$Method == "Fisher's Exact Test",c("term","term.y","pval","rank")] %>%
    as.data.frame()
  to_file_hyper <- to_file_hyper[order(to_file_hyper$pval),]
  colnames(to_file_hyper) <- c("Annotation Id","Term","P value","Rank")
  
  to_file_noncentral <- result_file[result_file$Method != "Fisher's Exact Test",c("term","term.y","pval","rank")] %>%
    as.data.frame()
  to_file_noncentral <- to_file_noncentral[order(to_file_noncentral$pval),]
  colnames(to_file_noncentral) <- c("Annotation Id","Term","P value","Rank")
  
  
  options(java.parameters = "-Xmx8000m")
  message(glue("Writing table for {annotation}"))
  invisible(gc())
  if (annotation == "KEGG"){
    invisible(file.remove("Tables/Supplementary Table 6.xlsx"))
    write.xlsx(to_file_hyper, file = "Tables/Supplementary Table 6.xlsx",
               append = FALSE,
               sheetName = paste0(annotation, " Fisher's Exact Test"),
               row.names = F)
    write.xlsx(to_file_noncentral, file = "Tables/Supplementary Table 6.xlsx",
               append = T,
               sheetName = paste0(annotation, " Wallenius' Test"),
               row.names = F)
  } else{
    write.xlsx(to_file_hyper, file = "Tables/Supplementary Table 6.xlsx",
               append = T,
               sheetName = paste0(annotation, " Fisher's Exact Test"),
               row.names = F)
    write.xlsx(to_file_noncentral, file = "Tables/Supplementary Table 6.xlsx",
               append = T,
               sheetName = paste0(annotation, " Wallenius' Test"),
               row.names = F)
  }
  
  return(result_file)
}))


sle_results <- rbindlist(lapply(annotations, function(annotation){
  result_file <- paste0("case_of_use/Cancer_",gsub(" ","_",annotation),".tsv")
  result_file <- vroom(result_file) %>%
    as.data.frame()
  
  result_file <- rbindlist(lapply(c("hyper","noncentral"), function(method){
    columns <- grepl(method,colnames(result_file))
    columns <- colnames(result_file)[columns]
    content <- result_file[,c("term",columns)]
    rank_vals <- content[,str_detect(colnames(content),"rank")]
    pval <- content[,method]
    terms <- content$term
    content <- data.frame(term = terms, Method = method, rank = rank_vals, pval = pval)
  }))
  result_file$annotation <- annotation
  result_file$Method <- ifelse(result_file$Method == "hyper",
                               "Fisher's Exact Test","Wallenius' Test")
  result_file <- left_join(result_file, ann_info, by = c("term" = "annotation_id"))
  to_file_hyper <- result_file[result_file$Method == "Fisher's Exact Test",c("term","term.y","pval","rank")] %>%
    as.data.frame()
  to_file_hyper <- to_file_hyper[order(to_file_hyper$pval),]
  colnames(to_file_hyper) <- c("Annotation Id","Term","P value","Rank")
  
  to_file_noncentral <- result_file[result_file$Method != "Fisher's Exact Test",c("term","term.y","pval","rank")] %>%
    as.data.frame()
  to_file_noncentral <- to_file_noncentral[order(to_file_noncentral$pval),]
  colnames(to_file_noncentral) <- c("Annotation Id","Term","P value","Rank")
  
  
  options(java.parameters = "-Xmx8000m")
  message(glue("Writing table for {annotation}"))
  invisible(gc())
  if (annotation == "KEGG"){
    invisible(file.remove("Tables/Supplementary Table 7.xlsx"))
    write.xlsx(to_file_hyper, file = "Tables/Supplementary Table 7.xlsx",
               append = FALSE,
               sheetName = paste0(annotation, " Fisher's Exact Test"),
               row.names = F)
    write.xlsx(to_file_noncentral, file = "Tables/Supplementary Table 7.xlsx",
               append = T,
               sheetName = paste0(annotation, " Wallenius' Test"),
               row.names = F)
  } else{
    write.xlsx(to_file_hyper, file = "Tables/Supplementary Table 7.xlsx",
               append = T,
               sheetName = paste0(annotation, " Fisher's Exact Test"),
               row.names = F)
    write.xlsx(to_file_noncentral, file = "Tables/Supplementary Table 7.xlsx",
               append = T,
               sheetName = paste0(annotation, " Wallenius' Test"),
               row.names = F)
  }
  
  return(result_file)
}))
