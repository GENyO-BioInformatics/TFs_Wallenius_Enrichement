library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(colorspace)
library(vroom)
library(parallel)
library(ggpubr)
library(ggdist)
library(rjson)
library(sdamr)
library(ggridges)
library(ggnewscale)
library(ggpp)
library(see)
library(tidyverse)

library(ggparl)

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

annotations = c("KEGG","Reactome","GO_BP","WikiPathways")

#https://medialab.github.io/iwanthue/
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

ann_info <- vroom("data/annotation_info_table.tsv")

### Plot Gene-Term Pairs (%). Figure 2.tiff

collectriTFsGRN <- read.delim('data/collectri_raw.tsv',header = T,sep = '\t')

# trrust <- vroom("trrust_rawdata.human.tsv",col_names = F)[,1:2]
# colnames(trrust) <- c("TF","target")

pairs_term_gene <- do.call("rbind",lapply(annotations, function(annotation){
  data_content <- vroom(paste0("data/",annotation,".tsv")) %>%
    filter(organism == 9606)
  
  ## Number Gene-Term Pairs
  n <- nrow(unique(data_content))
  
  tfs <- unique(dorothea_info$tf)
  
  only_tfs <- data_content[data_content$symbol %in% tfs,]
  n_genes <- length(unique(only_tfs$symbol))
  only_tfs <- round((nrow(unique(only_tfs)) / n) * 100,2)
  only_tfs <- data.frame("term_gene" = only_tfs, "n_genes" = n_genes)
  only_tfs$confidence <- "DoRothEA TFs"
  only_tfs$type = "TFs"
  
  confidence_level <- do.call("rbind",mclapply(c("A","B","C","D","E"), function(confidence){
    dorothea_filt <- dorothea_info[dorothea_info$confidence <= confidence,]
    target_genes <- unique(dorothea_filt$target)
    target_genes <- data_content[data_content$symbol %in% target_genes,]
    n_genes <- length(unique(target_genes$symbol))
    target_genes <- round((nrow(unique(target_genes))/n)*100,2)
    target_genes <- data.frame("term_gene" = target_genes, "n_genes" = n_genes, "confidence" = paste0("DoRothEA ",confidence), "type" = "Target Genes")
  }))
  confidence_level <- rbind(confidence_level,only_tfs)
  
  tfs <- unique(trrust$TF)
  
  only_tfs <- data_content[data_content$symbol %in% tfs,]
  n_genes <- length(unique(only_tfs$symbol))
  only_tfs <- round((nrow(unique(only_tfs)) / n) * 100,2)
  only_tfs <- data.frame("term_gene" = only_tfs, "n_genes" = n_genes)
  only_tfs$confidence <- "TRRUST v2 TFs"
  only_tfs$type = "TFs"
  
  confidence_level <- rbind(confidence_level,only_tfs)
  
  target_genes <- unique(trrust$target)
  target_genes <- data_content[data_content$symbol %in% target_genes,]
  n_genes <- length(unique(target_genes$symbol))
  target_genes <- round((nrow(unique(target_genes))/n)*100,2)
  
  confidence_level <- rbind(confidence_level,
                            data.frame("term_gene" = target_genes, 
                                       "n_genes" = n_genes, 
                                       "confidence" = "TRRUST v2",
                                       "type" = "Target Genes"))
  
  
  confidence_level$annotation <- gsub("_"," ",annotation)
  return(confidence_level)
}))

pairs_term_gene_file = pairs_term_gene[,c("annotation","confidence","term_gene","n_genes","type")]
colnames(pairs_term_gene_file) <- c("Annotation","TFs / Target confidence", "% coverage annotation", "Number of genes","Type")

library(xlsx)

write.xlsx(pairs_term_gene_file, file = "Supplementary Table 1.xlsx",
           append = FALSE,
           row.names = F)

annotations_names <- gsub("_"," ",annotations)

pairs_term_gene$annotation <- factor(pairs_term_gene$annotation, levels = annotations_names)

pairs_term_gene <- pairs_term_gene[order(pairs_term_gene$term_gene),]
pairs_term_gene$confidence <- factor(pairs_term_gene$confidence, levels = unique(pairs_term_gene$confidence))

pairs_term_gene$annotation <- factor(pairs_term_gene$annotation, levels = gsub("_"," ",annotations))

pairs_term_gene$type <- factor(pairs_term_gene$type, levels = c("TFs","Target Genes"))

gg <- ggplot(pairs_term_gene, aes(y = term_gene, x = confidence, colour = annotation))+
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

ggsave("Figure 2.tiff", plot = gg,units = "cm",height = 8, width = 8,dpi = 500)


### Raincloud + BarPlot Figure 3

gradient_colors <- lapply(names(colors_blind), function(annotation){
  base_color <- colors_blind[[annotation]]
  dark_color <- darken(base_color, 0.2)
  light_color <- lighten(base_color,0.6)
  colors_vec <- colorRampPalette(c(light_color,dark_color))(10)
})

names(gradient_colors) <- names(colors_blind)

dorothea_info <- vroom("data/dorothea.tsv") %>%
  filter(org == 9606, confidence <= "C")

tfs <- unique(dorothea_info$tf)

count_genes <- do.call("rbind",lapply(annotations, function(annotation){
  data_content <- vroom(paste0("data/",annotation,".tsv")) %>%
    filter(organism == 9606)
  target_genes <- unique(dorothea_info$target)
  target_genes <- data_content[data_content$symbol %in% target_genes,]
  target_genes <- unique(target_genes$symbol)
  dorothea_filt <- dorothea_info[dorothea_info$target %in% target_genes,] %>%
    group_by(target) %>%
    summarise(tfs = n()) %>%
    mutate(annotation = gsub("_"," ",annotation))
  return(dorothea_filt)
}))

count_genes_top <- count_genes %>%
  group_by(annotation) %>%
  top_n(10,wt=tfs)

count_genes_top <- unique(count_genes_top$target)

count_genes$annotation <- factor(count_genes$annotation, levels = rev(unique(count_genes$annotation)))

gg1 <- ggplot(count_genes, aes(annotation,y = tfs,fill=annotation))+
  geom_flat_violin(position = position_nudge(x = .05, y = 0), 
                   alpha = .8, adjust = 2) +
  geom_point(aes(y = tfs, color = annotation), position = position_jitternudge(nudge.x = -0.3,jitter.width = 0.2,seed = 123), size = .5, alpha = 0.8,
             adjust = 2) +
  geom_boxplot(width = .3, guides = FALSE, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = -.3, y = 0))+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "none",
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 6))+
  scale_color_manual(values = colors_blind) +
  scale_fill_manual(values = colors_blind) +
  ylab("TFs that regulate each gene")+
  geom_text_repel(data = count_genes[count_genes$target %in% count_genes_top,],aes(annotation,tfs,label=target),
                  position = position_jitternudge(nudge.x = -0.3,jitter.width = 0.1,seed = 123),size = 1, min.segment.length = 0.3)+
  coord_flip()

ann_info <- vroom("data/annotation_info_table.tsv")

ann_table <- do.call("rbind",lapply(annotations, function(annotation){
  count_genes <- count_genes[count_genes$annotation == gsub("_"," ",annotation),] %>%
    arrange(desc(tfs)) %>%
    top_n(20, wt = tfs)
  ann_table <- paste0("data/",annotation,".tsv") %>%
    vroom() %>%
    filter(organism == 9606, symbol %in% count_genes$target)
  ann_table_show <- as.data.frame(table(ann_table$annotation_id)) %>%
    arrange(desc(Freq)) %>%
    top_n(10,wt = Freq)
  ann_table_show$annotation_id <- as.character(ann_table_show$Var1)
  ann_table_show$annotation <- gsub("_"," ",annotation)
  ann_table_show <- merge(ann_table_show,ann_info,by = "annotation_id") %>%
    select(annotation_id,annotation,term,Freq)
  return(ann_table_show)
}))

ann_table <- ann_table[order(ann_table$Freq,decreasing = T),]
ann_table$annotation_id <- factor(ann_table$annotation_id,levels = unique(ann_table$annotation_id))

gg2 <- ggplot(ann_table, aes(x = Freq, y = annotation_id, fill = annotation, label = term)) +
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
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(0.3,"cm"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  labs(fill = "Annotation")+
  xlab("Target genes related to each term")

gg <- plot_grid(gg1,gg2,align = "h")

ggsave("Figure 3.tiff", plot = gg,units = "cm",height = 12, width = 20,dpi = 500)

library(parallel)

content <- do.call("rbind", lapply(annotations, function(annotation){
  result_files = list.files("null_simulation_tfs_rank/",pattern = paste0(annotation,"_"))
  nsize = unlist(strsplit(result_files,".tsv"))
  nsize = as.numeric(sapply(strsplit(nsize,paste0(annotation,"_")), `[`, 2))
  nsize = nsize[order(nsize)]
  
  content <- do.call("rbind",mclapply(nsize, function(n){
    result_file = paste0("null_simulation_tfs_rank/",annotation,"_",n,".tsv")
    content = vroom::vroom(result_file) %>%
      group_by(term) %>%
      summarise(mean_rank_hypergeom = mean(rank_hyper),
                proportion = sum(pval_hyper < 0.05) / 1000)
    content$n = n
    return(content)
  }))
  
  content$n <- factor(content$n, levels = unique(content$n))
  content$term <- factor(content$term, levels = unique(content$term))
  content$annotation = gsub("_"," ",annotation)
  return(content)
}))


gg1_1 <- ggplot(content, aes(x = n, y = mean_rank_hypergeom, color = annotation))+
  # geom_boxjitter(position = position_dodge(width = 1), jitter.size = 0.05, outlier.shape = NA, lwd = 0.2)+
  geom_boxplot(outlier.size = 0.2, lwd = 0.2)+
  scale_color_manual(values = colors_blind)+
  scale_y_log10()+
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
  ylab("Mean Ranking")+
  coord_flip()

gg1_2 <- ggplot(content, aes(x = n, y = proportion, color = annotation))+
  # geom_boxjitter(position = position_dodge(width = 1), jitter.size = 0.05, outlier.shape = NA, lwd = 0.2)+
  geom_boxplot(outlier.size = 0.2, lwd = 0.2)+
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
  ylab("Proportion of times")+
  coord_flip()


content_top = content %>%
  group_by(annotation) %>%
  top_n(-10,wt = mean_rank_hypergeom)

content_top = content[content$term %in% content_top$term,]
content_top = content_top[order(content_top$n),]
content_top$n <- paste0("N TFs: ",content_top$n)
content_top$n <- factor(content_top$n, levels = unique(content_top$n))

gg2_1 <- ggplot(content_top,aes(x = mean_rank_hypergeom, y = term, fill = annotation))+
  facet_wrap(~n, nrow = 1)+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = colors_blind)+
  scale_x_log10()+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title.x  = element_text(size = 5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4))+
  guides(fill = guide_legend(title = "Annotation"))+
  xlab("Mean Ranking")

gg_legend <- get_legend(gg2_1)

gg2_1 <- gg2_1 + theme(legend.position = "none")

gg_1 <- plot_grid(gg1_1, gg2_1, nrow = 2, rel_heights = c(0.8,1.2))

gg_1 <- plot_grid(gg_1, gg_legend, rel_widths = c(0.9,0.2))

ggsave("Supplementary Figure 1.tiff", plot = gg_1,units = "cm",height = 18, width = 20,dpi = 500)  

content_top = content %>%
  group_by(annotation) %>%
  top_n(10,wt = proportion)

content_top = content[content$term %in% content_top$term,]
content_top = content_top[order(content_top$n),]
content_top$n <- paste0("N TFs: ",content_top$n)
content_top$n <- factor(content_top$n, levels = unique(content_top$n))

gg2_2 <- ggplot(content_top,aes(x = proportion, y = term, fill = annotation))+
  facet_wrap(~n, nrow = 1)+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = colors_blind)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 4,face = "bold"),
        axis.title.x  = element_text(size = 5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4))+
  guides(fill = guide_legend(title = "Annotation"))+
  xlab("Mean Ranking")

gg_legend <- get_legend(gg2_2)

gg2_2 <- gg2_2 + theme(legend.position = "none")

gg_2 <- plot_grid(gg1_2, gg2_2, nrow = 2, rel_heights = c(0.8,1.2))

gg_2 <- plot_grid(gg_2, gg_legend, rel_widths = c(0.9,0.2))

ggsave("Supplementary Figure 1 v2.tiff", plot = gg_2,units = "cm",height = 18, width = 20,dpi = 500)  

### Plot Null Sim Target Genes results ###
annotation = "KEGG"

content <- lapply(annotations, function(annotation){
  nsize = c(2,3,4,5,10,20,30)
  
  content <- mclapply(nsize, function(n){
    result_file = paste0("null_simulation_target_genes_rank/",annotation,"_",n,".tsv")
    full_content <- vroom(result_file)
    
    content = full_content %>%
      group_by(term) %>%
      summarise(median_rank_hypergeom = median(rank_hyper),
                median_rank_nfisher = median(rank_nfisher),
                proportion_hypergeom = sum(hyper < 0.05) / 1000,
                proportion_nfisher = sum(nfisher < 0.05) / 1000)
    
    content <- do.call("rbind",lapply(c("hypergeom","nfisher"), function(method){
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
                             "Fisher","Fisher non-central")
    
    full_content <- do.call("rbind",lapply(c("hyper","nfisher"), function(method){
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
                                  "Fisher","Fisher non-central")
    content <- list(content,full_content)
    return(content)
  })
  
  full_content <- do.call("rbind",sapply(content,`[`,2)) 
  content <- do.call("rbind",sapply(content,`[`,1)) 
  
  content$n <- factor(content$n, levels = unique(content$n))
  content$annotation = gsub("_"," ",annotation)
  full_content$n <- factor(full_content$n, levels = unique(full_content$n))
  full_content$annotation = gsub("_"," ",annotation)
  content = list(content, full_content)
  return(content)
})

full_content <- do.call("rbind",sapply(content,`[`,2))
content <- do.call("rbind",sapply(content,`[`,1))

full_content <- full_content[full_content$term %in% ann_info$annotation_id,]
content <- content[content$term %in% ann_info$annotation_id,]

color_vals <- c()
for (annotation in annotations){
  base_color <- colors_blind[names(colors_blind) %in% gsub("_"," ",annotation)]
  dark_color <- darken(base_color, 0.3)
  light_color <- lighten(base_color,0.3)
  color_vals_insert <- c(dark_color, light_color)
  names(color_vals_insert) <- c(paste0(gsub("_"," ",annotation)," Fisher"),
                                paste0(gsub("_"," ",annotation)," Fisher non-central"))
  color_vals <- c(color_vals, color_vals_insert)
}

content$intersect <- paste0(gsub("_"," ",content$annotation)," ",content$Method)
content$intersect <- factor(content$intersect, levels = unique(content$intersect))

full_content$intersect <- paste0(gsub("_"," ",full_content$annotation)," ",full_content$Method)
full_content$intersect <- factor(full_content$intersect, levels = unique(full_content$intersect))

content$annotation <- factor(content$annotation, levels = gsub("_"," ",annotations))

gg_1_1 <- ggplot(content, aes(x = n, y = median_rank, color = intersect))+
  facet_wrap(~annotation, scales = "free")+
  geom_boxplot(outlier.size = 0.2, lwd = 0.2)+
  scale_color_manual(values = color_vals)+
  scale_y_log10()+
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
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5))+
  guides(color = guide_legend(title = "Annotation - Method"))+
  xlab("Number of TFs")+
  ylab("Median Ranking")


top_terms <- content %>%
  filter(n == 30, Method == "Fisher") %>%
  group_by(annotation) %>%
  top_n(-3,wt = median_rank)


sub_table <- full_content[full_content$term %in% top_terms$term,]
sub_table <- merge(sub_table,ann_info, by.x = "term", by.y = "annotation_id")
sub_table$annotation <- factor(sub_table$annotation, levels = gsub("_"," ",annotations))
sub_table <- sub_table[order(sub_table$annotation),]
sub_table$term_ann <- paste0(sub_table$annotation," - ",sub_table$term.y)
sub_table$term_ann <- factor(sub_table$term_ann, levels = unique(sub_table$term_ann))


library(ggpubr)
gg_2_1 <- ggplot(sub_table, aes(x = n, y = rank, color = intersect))+
  facet_wrap(~ term_ann,nrow = 4,scales = "free") +
  geom_boxplot(outlier.size = 0.3, lwd = 0.3)+
  stat_compare_means(method = "wilcox", label = "p.format", size = 1)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 3,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.position = "none")+
  scale_color_manual(values = color_vals)+
  xlab("Number of TFs")+
  ylab("Median Ranking")

gg_legend <- get_legend(gg_1_1)
gg_1_1 <- gg_1_1 + theme(legend.position = "none")
gg <- plot_grid(gg_1_1, gg_2_1, nrow = 2, rel_heights = c(0.8,1.2))
gg <- plot_grid(gg, gg_legend, rel_widths = c(0.9,0.2))
ggsave("Figure 4.tiff", plot = gg,units = "cm",height = 18, width = 20,dpi = 500)

############# Proportions ###################

gg_1_2 <- ggplot(content, aes(x = n, y = proportion, color = intersect))+
  facet_wrap(~annotation, scales = "free", nrow = 1)+
  geom_boxplot(outlier.size = 0.2, lwd = 0.2)+
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
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5))+
  guides(color = guide_legend(title = "Annotation - Method"))+
  xlab("Number of TFs")+
  ylab("Proportion of times significant (p< 0.05)")


top_terms <- content %>%
  filter(n == 30, Method == "Fisher") %>%
  arrange(desc(proportion)) %>%
  group_by(annotation) %>%
  do(head(.,3))


sub_table <- content[content$term %in% top_terms$term,]
sub_table <- merge(sub_table,ann_info, by.x = "term", by.y = "annotation_id")
sub_table$annotation <- factor(sub_table$annotation, levels = gsub("_"," ",annotations))
sub_table <- sub_table[order(sub_table$annotation),]
sub_table$term_ann <- paste0(sub_table$annotation," - ",sub_table$term.y)
sub_table$term_ann <- factor(sub_table$term_ann, levels = unique(sub_table$term_ann))


library(ggpubr)
gg_2_2 <- ggplot(sub_table, aes(x = n, y = proportion, color = intersect))+
  facet_wrap(~ term_ann, scales = "free", nrow = 4) +
  geom_point(size = 0.5)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 4),
        title = element_text(face = "bold"),
        strip.text = element_text(size = 3,face = "bold"),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.position = "none")+
  scale_color_manual(values = color_vals)+
  xlab("Number of TFs")+
  ylab("Proportion of times significant (p< 0.05)")

gg_legend <- get_legend(gg_1_2)
gg_1_2 <- gg_1_2 + theme(legend.position = "none")
gg <- plot_grid(gg_1_2, gg_2_2, nrow = 2)
gg <- plot_grid(gg, gg_legend, rel_widths = c(0.9,0.2))
ggsave("Figure 4 v2.tiff", plot = gg,units = "cm",height = 18, width = 20,dpi = 500)



content <- lapply(annotations,function(annotation){
  result_files = list.files("null_simulation_target_genes_rank/",pattern = paste0(annotation,"_"))
  nsize = unlist(strsplit(result_files,".tsv"))
  nsize = as.numeric(sapply(strsplit(nsize,paste0(annotation,"_")), `[`, 2))
  nsize = nsize[order(nsize)]
  
  content <- do.call("rbind",mclapply(nsize, function(n){
    result_file = paste0("null_simulation_target_genes_rank/",annotation,"_",n,".tsv")
    content = vroom::vroom(result_file) %>%
      group_by(term) %>%
      summarise(mean_rank_hypergeom = mean(rank_hyper),
                mean_rank_wallenius = mean(rank_wall)) %>%
      pivot_longer(cols = c("mean_rank_hypergeom",
                            "mean_rank_wallenius"),
                   names_to = "method",
                   values_to = "value")
    
    content <- do.call("rbind",lapply(c("hypergeom","wallenius"), function(method){
      rows <- grepl(method,content$method)
      content <- content[rows,]
      mean_rank <- content[str_detect(content$method,"mean"),]$value
      terms <- unique(content$term)
      content <- data.frame(term = terms, Method = method, mean_rank = mean_rank)
    }))
    
    content$n = n
    
    content$Method <- ifelse(content$Method == "hypergeom",
                             "Fisher","Wallenius")
    return(content)
  }))
  
  content$n <- paste0("N TFs: ",content$n)
  content$n <- factor(content$n, levels = unique(content$n))
  order_terms <- content[content$Method == "Fisher",]
  order_terms <- order_terms[order(order_terms$mean_rank),"term"]
  content$term <- factor(content$term, levels = unique(order_terms))
  
  base_color <- colors_blind[names(colors_blind) %in% gsub("_"," ",annotation)]
  dark_color <- darken(base_color, 0.3)
  light_color <- lighten(base_color,0.3)
  color_vals <- c("Fisher" = dark_color, "Wallenius" = light_color)
  
  dark_color <- base_color
  names(dark_color) <- NULL
  light_color <- lighten(base_color,0.4)
  color_vals_points <- c("Fisher" = dark_color, "Wallenius" = light_color)
  
  
  ggplot(content, aes(y = mean_rank, x = n, color = Method))+
    geom_jitter(position = position_dodge())
  
  gg <- ggplot(content, aes(y = mean_rank, x = term, group = Method,
                            color = Method))+
    facet_wrap(~n, nrow = 4)+
    geom_smooth(size = 0.5, method = "loess",se = F)+
    scale_color_manual(values = color_vals)+
    theme(panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          legend.key = element_rect(fill = "white"),
          title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 6),
          strip.text = element_text(size = 4,face = "bold"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          axis.title = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.key.size  = unit(0.3, 'cm'),
          plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
          legend.margin=margin(0,8,0,0),
          legend.box.margin=margin(-10,0,-10,-10))+
    ylab("Terms")+
    xlab("Mean Ranking")+
    ggtitle(gsub("_"," ",annotation))
  # ggsave(paste0(annotation,".tiff"), plot = gg,units = "cm",height = 10, width = 15,dpi = 500)
  return(gg)
})

big_plot <- plot_grid(plotlist = content,align = "hv")
ggsave("Figure 4.tiff", plot = big_plot,units = "cm",height = 12, width = 20,dpi = 500)


content <- lapply(annotations,function(annotation){
  result_files = list.files("null_simulation_target_genes_rank/",pattern = paste0(annotation,"_"))
  nsize = unlist(strsplit(result_files,".tsv"))
  nsize = as.numeric(sapply(strsplit(nsize,paste0(annotation,"_")), `[`, 2))
  nsize = nsize[order(nsize)]
  
  content <- do.call("rbind",mclapply(nsize, function(n){
    result_file = paste0("null_simulation_target_genes_rank/",annotation,"_",n,".tsv")
    content = vroom::vroom(result_file) %>%
      group_by(term) %>%
      summarise(mean_rank_hypergeom = mean(rank_hyper),
                mean_rank_wallenius = mean(rank_wall)) %>%
      pivot_longer(cols = c("mean_rank_hypergeom",
                            "mean_rank_wallenius"),
                   names_to = "method",
                   values_to = "value")
    
    content <- do.call("rbind",lapply(c("hypergeom","wallenius"), function(method){
      rows <- grepl(method,content$method)
      content <- content[rows,]
      mean_rank <- content[str_detect(content$method,"mean"),]$value
      terms <- unique(content$term)
      content <- data.frame(term = terms, Method = method, mean_rank = mean_rank)
    }))
    
    content$n = n
    
    content$Method <- ifelse(content$Method == "hypergeom",
                             "Fisher","Wallenius")
    return(content)
  }))
  
  content$n <- paste0("N TFs: ",content$n)
  content$n <- factor(content$n, levels = unique(content$n))
  ann_table_sub <- ann_table[ann_table$annotation_id %in% content$term,]
  content <- content[content$term %in% ann_table_sub$annotation_id,]
  content$term <- factor(content$term, levels = unique(ann_table_sub$annotation_id))
  
  base_color <- colors_blind[names(colors_blind) %in% gsub("_"," ",annotation)]
  dark_color <- darken(base_color, 0.3)
  light_color <- lighten(base_color,0.3)
  color_vals <- c("Fisher" = dark_color, "Wallenius" = light_color)
  
  
  gg <- ggplot(content, aes(x = mean_rank, y = term, fill = Method))+
    facet_wrap(~n, nrow = 4)+
    geom_bar(stat = "identity", position = "dodge")+
    scale_fill_manual(values = color_vals)+
    theme(panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 4),
          legend.key = element_rect(fill = "white"),
          title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 6),
          strip.text = element_text(size = 4,face = "bold"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          axis.title = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.key.size  = unit(0.3, 'cm'),
          plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
          legend.margin=margin(0,8,0,0),
          legend.box.margin=margin(-10,0,-10,-10))+
    ylab("Terms")+
    xlab("Mean Ranking")+
    ggtitle(gsub("_"," ",annotation))
  
  
  
  # gg <- ggplot(content, aes(x = value, y = term,
  #                           color = Method))+
  #   facet_wrap(~n, nrow = 4)+
  #   geom_boxplot(outlier.shape = NA) +
  #   scale_color_manual(values = color_vals)+
  #   theme(panel.background = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         strip.background = element_rect(fill = "white"),
  #         axis.text.x = element_text(size = 4),
  #         legend.key = element_rect(fill = "white"),
  #         title = element_text(face = "bold"),
  #         plot.title = element_text(hjust = 0.5,size = 6),
  #         strip.text = element_text(size = 4,face = "bold"),
  #         legend.text = element_text(size = 4),
  #         legend.title = element_text(size = 5),
  #         axis.title = element_text(size = 5),
  #         axis.text.y = element_text(size = 4),
  #         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         legend.key.size  = unit(0.3, 'cm'),
  #         plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
  #         legend.margin=margin(0,8,0,0),
  #         legend.box.margin=margin(-10,0,-10,-10))+
  #   ylab("Terms")+
  #   xlab("Ranking")+
  #   ggtitle(gsub("_"," ",annotation))
  ggsave(paste0(annotation,".tiff"), plot = gg,units = "cm",height = 10, width = 15,dpi = 500)
  return(gg)
})

big_plot <- plot_grid(plotlist = content,align = "hv")
ggsave("Figure 5.tiff", plot = big_plot,units = "cm",height = 18, width = 15,dpi = 500)

# 
# 
# 
# 
# library(network)
# library(sna)
# library(ggplot2)
# library(igraph)
# library(GGally)
# 
# genes <- c("MYC","CCND1","CDKN1A")
# 
# dorothea_info <- vroom("data/dorothea.tsv") %>%
#   filter(org == 9606, confidence <= "C", target %in% genes) %>%
#   select(tf,target) %>%
#   mutate(match = 1)
# 
# net <- get.adjacency(
#   graph_from_data_frame(dorothea_info, directed = FALSE),
#   attr = "match",
#   sparse = FALSE
# )
# 
# bip = network(net)
# 
# x = network.vertex.names(bip)
# x = ifelse(x %in% genes,"Target","TF")
# bip %v% "Type" = x
# 
# # detect and color the mode
# network_gg <- ggnet2(bip, color = "Type",label = TRUE,
#                      palette = c("Target" = "gold", "TF" = "grey"),
#                      size = 3,label.size = 1,legend.size = 3)
# 
# png("Supplementary Figure 2.png",width = 2000,height = 1000, res = 300)          
# gg <- ggplot(count_genes, aes(annotation,y = json_cont,fill=annotation))+
#   geom_flat_violin(position = position_nudge(x = .05, y = 0), alpha = .8, adjust = 2)+
#   geom_point(aes(y = json_cont, color = annotation), position = position_jitternudge(nudge.x = -0.3,jitter.width = 0.2,seed = 123), size = .5, alpha = 0.8,
#              adjust = 2)+
#   geom_boxplot(width = .3, guides = FALSE, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = -.3, y = 0))+
#   theme(panel.background = element_blank(),
#         strip.background = element_rect(fill = "white"),
#         axis.title.y = element_blank(),
#         legend.key = element_rect(fill = "white"),
#         legend.position = "none",
#         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.title.x = element_text(size = 6),
#         axis.text = element_text(size = 5))+
#   scale_color_manual(values = colors_blind) +
#   scale_fill_manual(values = colors_blind) +
#   ylab("TFs that regulate each gene")+
#   geom_text_repel(data = count_genes[count_genes$gene %in% count_genes_top,],aes(annotation,json_cont,label=gene),
#                   position = position_jitternudge(nudge.x = -0.3,jitter.width = 0.1,seed = 123),size = 1, min.segment.length = 0.3)+
#   coord_flip()
# gg <- plot_grid(gg,network_gg)
# plot(gg)
# invisible(dev.off())
# 
# count_genes_top <- count_genes %>%
#   group_by(annotation) %>%
#   top_n(20,wt=json_cont)
# 
# 
# 
# ann_table <- ann_table[,c("annotation","annotation_id","term","Freq")]
# colnames(ann_table) <- c("annotation","Annotation ID","Term","Count")
# 
# library(xlsx)

# for (annotation in annotations){
#   data_sub <- ann_table[ann_table$annotation == gsub("_"," ",annotation),c("Annotation ID","Term","Count")]
#   if (!(file.exists("supp_table_3.xlsx"))){
#     write.xlsx(data_sub, file = "supp_table_3.xlsx",
#                sheetName = annotation, append = FALSE,
#                row.names = F)
#   } else{
#     write.xlsx(data_sub, file = "supp_table_3.xlsx",
#                sheetName = annotation, append = TRUE,
#                row.names = F)
#   }
# }


## Plot null simulations using TFs


# content <- lapply(annotations,function(annotation){
#   result_files = list.files("null_simulation_tfs",pattern = annotation)
#   nsize = unlist(strsplit(result_files,".tsv"))
#   nsize = as.numeric(sapply(strsplit(nsize,paste0(annotation,"_")), `[`, 2))
#   nsize = nsize[order(nsize)]
#   
#   content <- do.call("rbind",lapply(nsize, function(n){
#     result_file = paste0("null_simulation_tfs/",annotation,"_",n,".tsv")
#     content = vroom(result_file) %>%
#       pivot_longer(cols = c("hypergeom"),names_to = "method",values_to = "proportion")
#     content$n = n
#     content$method <- ifelse(content$method == "hypergeom",
#                              "Fisher","Wallenius")
#     return(content)
#   }))
#   
#   content_table <- content[,c("term","method","proportion","n")]
#   colnames(content_table) <- c("Term","Method","Proportion of times significant (p < 0.05)","Number of TFs")
#   
#   content_table <- merge(content_table,ann_info,by.x = "Term",by.y = "annotation_id")
#   colnames(content_table)[1] <- "Annotation ID"
#   colnames(content_table)[ncol(content_table)] <- "Term"
#   content_table <- content_table[,c("Method","Annotation ID","Term","Number of TFs","Proportion of times significant (p < 0.05)")]
#   content_table <- content_table[order(content_table$`Proportion of times significant (p < 0.05)`,decreasing = T),]
#   
#   write.xlsx(content_table, file = paste0(gsub("_"," ",annotation),"_Supplementary Table 2.xlsx"),
#              sheetName = gsub("_"," ",annotation), append = FALSE,
#              row.names = F)
#   
#   content$annotation = gsub("_"," ",annotation)
#   content$n <- paste0("Number TFs: ",content$n)
#   content$n <- factor(content$n, levels = unique(content$n))
#   
#   # gg <- ggplot(data = content, aes(y = proportion, x = method, fill = method)) +
#   #   facet_wrap(~n,nrow = 1)+
#   #   geom_flat_violin(alpha = .8,position = position_nudge(x = .13, y = 0),lwd = 0.2,width = 0.8) +
#   #   geom_point(aes(y = proportion, color = method), position = position_dodgenudge(x = -.1, y = 0,width = 0.2), size = .05, alpha = 0.8) +
#   #   geom_boxplot(width = .2, guides = FALSE, outlier.shape = NA, alpha = 0.5, lwd = 0.2,position = position_nudge(x = -.1, y = 0)) +
#   #   guides(fill = "none")+
#   #   guides(colour = guide_legend(override.aes = list(size=1)))+
#   #   scale_fill_manual(values = color_vals)+
#   #   scale_color_manual(values = color_vals)+
#   #   theme(panel.background = element_blank(),
#   #         axis.ticks.x = element_blank(),
#   #         strip.background = element_rect(fill = "white"),
#   #         axis.title.x = element_blank(),
#   #         axis.text.x = element_blank(),
#   #         legend.key = element_rect(fill = "white"),
#   #         title = element_text(face = "bold"),
#   #         plot.title = element_text(hjust = 0.5,size = 3),
#   #         strip.text = element_text(size = 2,face = "bold"),
#   #         legend.text = element_text(size = 3),
#   #         legend.title = element_text(size = 3),
#   #         axis.title = element_text(size = 3),
#   #         axis.text.y = element_text(size = 3),
#   #         axis.line.x = element_line(colour = 'black', size=0.2, linetype='solid'),
#   #         axis.line.y = element_line(colour = 'black', size=0.2, linetype='solid'),
#   #         legend.key.size  = unit(0.3, 'cm'),
#   #         plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
#   #         legend.margin=margin(0,8,0,0),
#   #         legend.box.margin=margin(-10,0,-10,-10))+
#   #   ylab("Proportion of significance (p < 0.05)")+
#   #   labs(color = "Method")+
#   #   ggtitle(annotation)
#   # ggsave(paste0(annotation,".tiff"), plot = gg,units = "cm",height = 5, width = 8,dpi = 500)
#   # 
#   return(content)
# })

### Plot bias by gene ###



### Plot Null Sim Target Genes results ###
annotation = "KEGG"
content <- lapply(annotations,function(annotation){
  result_files = list.files("null_simulation_target_genes_rank/",pattern = paste0(annotation,"_"))
  nsize = unlist(strsplit(result_files,".tsv"))
  nsize = as.numeric(sapply(strsplit(nsize,paste0(annotation,"_")), `[`, 2))
  nsize = nsize[order(nsize)]
  
  content <- do.call("rbind",lapply(nsize, function(n){
    result_file = paste0("null_simulation_target_genes_rank/",annotation,"_",n,".tsv")
    content = vroom(result_file) %>%
      group_by(term) %>%
      summarise(mean_rank_hypergeom = mean(rank_hyper),
                mean_rank_wallenius = mean(rank_wall),
                var_rank_hypergeom = var(rank_hyper),
                var_rank_wallenius = var(rank_wall)) %>%
      pivot_longer(cols = c("mean_rank_hypergeom",
                            "var_rank_hypergeom",
                            "mean_rank_wallenius",
                            "var_rank_wallenius"),
                   names_to = "method",
                   values_to = "value")
    
    content <- do.call("rbind",lapply(c("hypergeom","wallenius"), function(method){
      rows <- grepl(method,content$method)
      content <- content[rows,]
      mean_rank <- content[str_detect(content$method,"mean"),]$value
      var_rank <- content[str_detect(content$method,"var"),]$value
      terms <- unique(content$term)
      content <- data.frame(term = terms, Method = method, var_rank = var_rank, mean_rank = mean_rank)
    }))
    
    content$n = n
    
    content$Method <- ifelse(content$Method == "hypergeom",
                             "Fisher","Wallenius")
    return(content)
  }))
  
  content$n <- paste0("N TFs: ",content$n)
  content$n <- factor(content$n, levels = unique(content$n))
  # content[content$method == "Wallenius","mean_rank"] <- -content[content$method == "Wallenius","mean_rank"]
  # content <- content[content$term %in% unique(content$term)[1:30],]
  order_terms <- content[content$Method == "Fisher",]
  order_terms <- order_terms[order(order_terms$mean_rank),"term"]
  content$term <- factor(content$term, levels = unique(order_terms))
  
  base_color <- colors_blind[names(colors_blind) %in% gsub("_"," ",annotation)]
  dark_color <- darken(base_color, 0.6)
  light_color <- darken(base_color,0.2)
  color_vals <- c("Fisher" = dark_color, "Wallenius" = light_color)
  
  dark_color <- base_color
  names(dark_color) <- NULL
  light_color <- lighten(base_color,0.2)
  color_vals_points <- c("Fisher" = dark_color, "Wallenius" = light_color)
  
  gg <- ggplot(content, aes(y = mean_rank, x = term, colour = Method,
                            group=Method, fill = Method))+
    facet_wrap(~n, nrow = 4)+
    geom_point(data = content, aes(y = mean_rank, x = term, color = Method),
               size = 0.00000005, alpha = 0.8) +
    scale_colour_manual(values = color_vals_points)
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.key = element_rect(fill = "white"),
        title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,size = 6),
        strip.text = element_text(size = 4,face = "bold"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 4),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.key.size  = unit(0.3, 'cm'),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        legend.margin=margin(0,8,0,0),
        legend.box.margin=margin(-10,0,-10,-10))+
    xlab("Terms")+
    ylab("Mean Ranking")+
    ggtitle(gsub("_"," ",annotation))+   
    new_scale_color()+
    geom_smooth(size = 0.3, method = "loess",se = F)+
    scale_color_manual(values = color_vals)
  
  ggsave(paste0(annotation,".tiff"), plot = gg,units = "cm",height = 10, width = 15,dpi = 500)
  return(gg)
})

big_plot <- plot_grid(plotlist = content)
ggsave("Figure 2.tiff", plot = big_plot,units = "cm",height = 10, width = 15,dpi = 500)

content <- lapply(annotations,function(annotation){
  result_files = list.files("null_simulation_target_genes_rank/",pattern = paste0(annotation,"_"))
  nsize = unlist(strsplit(result_files,".tsv"))
  nsize = as.numeric(sapply(strsplit(nsize,paste0(annotation,"_")), `[`, 2))
  nsize = nsize[order(nsize)]
  content <- do.call("rbind",lapply(nsize, function(n){
    result_file = paste0("null_simulation_target_genes_rank/",annotation,"_",n,".tsv")
    content = vroom(result_file) %>%
      pivot_longer(cols = c("rank_hyper","rank_wall"),
                   names_to = "method",
                   values_to = "value")
    content$n = n
    content$Method <- ifelse(content$method == "rank_hyper",
                             "Fisher","Wallenius")
    return(content)
  }))
  
  content <- content[content$term %in% ann_table$`Annotation ID`,c("term","tfs","value","n","Method")]
  content <- merge(content, ann_table, by.x = "term", by.y = "Annotation ID")
  
  content <- content[order(content$n),]
  content$n <- paste0("N TFs: ",content$n)
  content$n <- factor(content$n, levels = unique(content$n))
  content <- content[order(content$value),]
  content$term <- factor(content$term, levels = unique(content$term))
  
  quantiles <- content %>%
    group_by(Method, term) %>%
    summarize(quant25 = quantile(value, probs = 0.75))
  
  max_val <- max(quantiles$quant25) + 10
  
  gg <- ggplot(content, aes(x = value, y = term))+
    geom_boxplot(aes(fill = Method), alpha = 0.5,outlier.shape = NA)+
    geom_text(data = unique(content[,c("term","Term")]),
              aes(x = max_val, y = term, label = Term),
              hjust = 0, size = 1.5)+
    scale_fill_manual(values = color_vals)+
    theme(panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.key = element_rect(fill = "white"),
          title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 6),
          strip.text = element_text(size = 4,face = "bold"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          axis.title = element_text(size = 5),
          axis.text.x = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.key.size  = unit(0.3, 'cm'),
          plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
          legend.margin=margin(0,8,0,0),
          legend.box.margin=margin(-10,0,-10,-10))+
    ylab("Terms")+
    xlab("Ranking")+
    ggtitle(gsub("_"," ",annotation))
  
  # ggsave(paste0(annotation,"_biased_terms.tiff"), plot = gg,units = "cm",height = 7, width = 15,dpi = 500)
  
  # 
  # 
  # 
  # content_table <- content[,c("term","method","proportion","n")]
  # colnames(content_table) <- c("Term","Method","Proportion of times significant (p < 0.05)","Number of TFs")
  # 
  # content_table <- merge(content_table,ann_info,by.x = "Term",by.y = "annotation_id")
  # colnames(content_table)[1] <- "Annotation ID"
  # colnames(content_table)[ncol(content_table)] <- "Term"
  # content_table <- content_table[,c("Method","Annotation ID","Term","Number of TFs","Proportion of times significant (p < 0.05)")]
  # 
  # hyper <- content_table[content_table$Method == "Fisher",-1]
  # hyper <- hyper[order(hyper$`Proportion of times significant (p < 0.05)`,decreasing = T),]
  # wall <- content_table[content_table$Method != "Fisher",-1]
  # wall <- wall[order(wall$`Proportion of times significant (p < 0.05)`,decreasing = T),]
  # 
  # 
  # # write.xlsx(hyper, file = paste0(gsub("_"," ",annotation),"_Supplementary Table 4.xlsx"),
  # #            sheetName = gsub("_"," ",annotation), append = FALSE,
  # #            row.names = F)
  # # write.xlsx(wall, file = paste0(gsub("_"," ",annotation),"_Supplementary Table 5.xlsx"),
  # #            sheetName = gsub("_"," ",annotation), append = FALSE,
  # #            row.names = F)
  # 
  # content$annotation = gsub("_"," ",annotation)
  # content$n <- paste0("N TFs: ",content$n)
  # content$n <- factor(content$n, levels = unique(content$n))
  # 
  # base_color <- colors_blind[names(colors_blind) %in% content$annotation]
  # dark_color <- darken(base_color, 0.2)
  # light_color <- lighten(base_color,0.6)
  # color_vals <- c("Fisher" = dark_color, "Wallenius" = light_color)
  # 
  # gg <- ggplot(data = content, aes(y = proportion, x = method, fill = method)) +
  #   facet_wrap(~n,nrow = 1)+
  #   geom_flat_violin(alpha = .8,position = position_nudge(x = .13, y = 0),lwd = 0.2,width = 0.8) +
  #   geom_point(aes(y = proportion, color = method), position = position_dodgenudge(x = -.1, y = 0,width = 0.2), size = .05, alpha = 0.8) +
  #   geom_boxplot(width = .2, guides = FALSE, outlier.shape = NA, alpha = 0.5, lwd = 0.2,position = position_nudge(x = -.1, y = 0)) +
  #   guides(fill = "none")+
  #   guides(colour = guide_legend(override.aes = list(size=2)))+
  #   scale_fill_manual(values = color_vals)+
  #   scale_color_manual(values = color_vals)+
  #   theme(panel.background = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         strip.background = element_rect(fill = "white"),
  #         axis.title.x = element_blank(),
  #         axis.text.x = element_blank(),
  #         legend.key = element_rect(fill = "white"),
  #         title = element_text(face = "bold"),
  #         plot.title = element_text(hjust = 0.5,size = 6),
  #         strip.text = element_text(size = 4,face = "bold"),
  #         legend.text = element_text(size = 4),
  #         legend.title = element_text(size = 5),
  #         axis.title = element_text(size = 5),
  #         axis.text.y = element_text(size = 4),
  #         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         legend.key.size  = unit(0.3, 'cm'),
  #         plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
  #         legend.margin=margin(0,8,0,0),
  #         legend.box.margin=margin(-10,0,-10,-10))+
  #   ylab("Proportion of significance (p < 0.05)")+
  #   labs(color = "Method")+
  #   ggtitle(annotation)
  # # ggsave(paste0(annotation,".tiff"), plot = gg,units = "cm",height = 5, width = 8,dpi = 500)
  # 
  return(gg)
})

big_plot <- plot_grid(plotlist = content)

ggsave("Figure 2_1.tiff", plot = big_plot,units = "cm",height = 15, width = 19,dpi = 500)

#### Use of Case SLE

ann_info <- vroom("data/annotation_info_table.tsv")
library(grid)

annotation <- "WikiPathways"
content <- lapply(annotations, function(annotation){
  base_color <- colors_blind[names(colors_blind) %in% gsub("_"," ",annotation)]
  dark_color <- darken(base_color, 0.2)
  light_color <- lighten(base_color,0.6)
  color_vals <- c("Fisher" = dark_color, "Wallenius" = light_color)
  result_file <- paste0("results/SLE_",annotation,".tsv")
  content <- vroom(result_file)
  
  
  
  
  content <- content[content$term %in% ann_table$`Annotation ID`,]
  pval_hyper <- content$pval_hyper
  pval_hyper <- data.frame(Term = content$term,Pval = pval_hyper, Method = "Fisher")
  pval_wall <- content$pval_wall
  pval_wall <- data.frame(Term = content$term,Pval = pval_wall, Method = "Wallenius")
  pval_comb <- rbind(pval_hyper,pval_wall)
  pval_comb <- merge(pval_comb,ann_info,by.x = "Term", by.y= "annotation_id") %>%
    arrange(Pval)
  colnames(pval_comb) <- c("ann_id","Pval","Method","Term")
  
  
  hyper <- pval_comb[pval_comb$Method == "Fisher",]
  wall <- pval_comb[pval_comb$Method == "Wallenius",]
  
  pval_comb$logPval <- -log(pval_comb$Pval)
  pval_comb <- pval_comb[order(pval_comb$logPval,decreasing = T),]
  pval_comb$Term <- factor(pval_comb$Term, levels = unique(pval_comb$Term))
  
  ggplot(pval_comb, aes(x = logPval, y = Term, fill = Method))+
    geom_bar(stat = "identity", position = position_dodge())+
    theme(panel.background = element_blank())+
    geom_vline(xintercept = -log(0.05))+
    geom_vline(xintercept = -log(0.01))+
    annotate("segment", x = 10, y = nrow(pval_comb)/2 - 5, xend = 3, yend = nrow(pval_comb)/2 - 5, 
             size = 0.3, linejoin = "mitre",
             arrow = arrow(type = "closed", length = unit(0.01, "npc")))+
    annotate("text", x = 10.5, y = nrow(pval_comb)/2 - 5, label = "pval 0.05",
             hjust = 0, size = 3)+
    annotate("segment", x = 10, y = nrow(pval_comb)/2 - 2, xend = 4.6, yend = nrow(pval_comb)/2 - 2, 
             size = 0.3, linejoin = "mitre",
             arrow = arrow(type = "closed", length = unit(0.01, "npc")))+
    annotate("text", x = 10.5, y = nrow(pval_comb)/2 - 2, label = "pval 0.01",
             hjust = 0, size = 3)+
    scale_fill_manual(values = color_vals)
  
  
  
  
  
  
  base_color <- colors_blind[names(colors_blind) %in% gsub("_"," ",annotation)]
  dark_color <- darken(base_color, 0.2)
  light_color <- lighten(base_color,0.4)
  color_vals <- c("Hypergeometric" = dark_color, "Wallenius" = light_color)
  
  (gg <- ggplot(pval_comb,aes(x = logPval,y = Term, label = Term, fill = Method))+
      scale_fill_manual(values = color_vals) +
      facet_wrap(~Method)+
      geom_bar(stat = "identity")+
      geom_text(vjust= 0.5, size = 1, hjust = 0, x = 0.5)+
      theme(axis.text.y = element_blank(),
            panel.background = element_blank(),
            strip.background = element_rect(fill = "white"),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.length.x = unit(0.05, "cm"),
            axis.ticks.x = element_line(size = 0.2),
            title = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5,size = 3),
            strip.text = element_text(size = 2,face = "bold"),
            legend.text = element_text(size = 2),
            legend.title = element_text(size = 3),
            axis.title = element_text(size = 3),
            legend.key.size  = unit(0.2, 'cm'),
            legend.key = element_rect(fill = "white"),
            plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
            legend.margin=margin(0,8,0,0),
            legend.box.margin=margin(-10,0,-10,-10),
            axis.text.x = element_text(size = 3),
            axis.line.x = element_line(colour = 'black', size=0.2, linetype='solid'),
            axis.line.y = element_line(colour = 'black', size=0.2, linetype='solid'))+
      ggtitle(gsub("_"," ",annotation))+
      xlab("- log (p-value)")+
      geom_vline(xintercept = 2.99, size = 0.2, colour = "darkred", linetype = 2)+
      geom_vline(xintercept = 4.6, size = 0.2, colour = "darkblue", linetype = 2))
  
  ggsave(paste0(annotation,".tiff"),units = "cm",height = 5, width = 8,dpi = 500)
})

content <- lapply(annotations, function(annotation){
  result_file <- paste0("results/Cancer_",annotation,".tsv")
  content <- vroom(result_file)
  pval_hyper <- content$pval_hyper
  pval_hyper <- data.frame(Term = content$term,Pval = pval_hyper, Method = "Hypergeometric")
  pval_wall <- content$pval_wall
  pval_wall <- data.frame(Term = content$term,Pval = pval_wall, Method = "Wallenius")
  pval_comb <- rbind(pval_hyper,pval_wall)
  pval_comb <- merge(pval_comb,ann_info,by.x = "Term", by.y= "annotation_id") %>%
    arrange(Pval)
  colnames(pval_comb) <- c("ann_id","Pval","Method","Term")
  
  sign_in_hyper <- pval_comb[pval_comb$Method == "Hypergeometric" & pval_comb$Pval < 0.05,"Term"]
  sign_in_wall <- pval_comb[pval_comb$Method == "Wallenius" & pval_comb$Pval < 0.05,"Term"]
  
  not_intersect <- sign_in_hyper[!(sign_in_hyper %in% sign_in_wall)]
  
  
  possible_terms <- na.omit(not_intersect[1:25])
  
  
  pval_comb <- pval_comb[pval_comb$Term %in% possible_terms,]
  
  terms <- pval_comb[pval_comb$Method == "Hypergeometric",]
  terms <- terms[order(terms$Pval),]
  pval_comb$Term <- factor(pval_comb$Term, levels = unique(terms$Term))
  pval_comb$logPval <- -log(pval_comb$Pval)
  
  base_color <- colors_blind[names(colors_blind) %in% gsub("_"," ",annotation)]
  dark_color <- darken(base_color, 0.2)
  light_color <- lighten(base_color,0.4)
  color_vals <- c("Hypergeometric" = dark_color, "Wallenius" = light_color)
  
  (gg <- ggplot(pval_comb,aes(x = logPval,y = Term, label = Term, fill = Method))+
      scale_fill_manual(values = color_vals) +
      facet_wrap(~Method)+
      geom_bar(stat = "identity")+
      geom_text(vjust= 0.5, size = 1, hjust = 0, x = 0.5)+
      theme(axis.text.y = element_blank(),
            panel.background = element_blank(),
            strip.background = element_rect(fill = "white"),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.length.x = unit(0.05, "cm"),
            axis.ticks.x = element_line(size = 0.2),
            title = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5,size = 3),
            strip.text = element_text(size = 2,face = "bold"),
            legend.text = element_text(size = 2),
            legend.title = element_text(size = 3),
            axis.title = element_text(size = 3),
            legend.key.size  = unit(0.2, 'cm'),
            legend.key = element_rect(fill = "white"),
            plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
            legend.margin=margin(0,8,0,0),
            legend.box.margin=margin(-10,0,-10,-10),
            axis.text.x = element_text(size = 3),
            axis.line.x = element_line(colour = 'black', size=0.2, linetype='solid'),
            axis.line.y = element_line(colour = 'black', size=0.2, linetype='solid'))+
      ggtitle(gsub("_"," ",annotation))+
      xlab("- log (p-value)")+
      geom_vline(xintercept = 2.99, size = 0.2, colour = "darkred", linetype = 2)+
      geom_vline(xintercept = 4.6, size = 0.2, colour = "darkblue", linetype = 2))
  
  ggsave(paste0(annotation,".tiff"),units = "cm",height = 5, width = 8,dpi = 500)
})
