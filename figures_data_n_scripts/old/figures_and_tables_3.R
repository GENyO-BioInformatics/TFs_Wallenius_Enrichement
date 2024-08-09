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
library(lemon)
library(gtable)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

### Prepare data

annotations = c("KEGG","Reactome","GO BP","WikiPathways")

colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "WikiPathways" = "#D974A0",
                  "Reactome" = "#f1b620")

ann_info <- vroom("data/annotation_info_table.tsv") %>%
  as.data.frame()

ann_tables <- rbindlist(mclapply(annotations, function(annotation){
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

regulons <- rbind(dorothea_info)

################################ Plot Collectri Annotations Coverage By Evidence 
################################ Figure 2.tiff

annotationsCoveragePerEvidence <- vroom("data/observations/annotationsCoveragePerEvidence.tsv")

annotationsCoveragePerEvidence$minEvidenceScore
annotationsCoveragePerEvidence$propor_TargetsAnnotationsInEvidence

annotationsCoveragePerEvidence$evidence

annotationsCoveragePerEvidence$evidence

color_vals2 <- colors_blind
names(color_vals2) <- annotations


color_vals[grep('Target',names(color_vals))]

byCuration_effortDF <- annotationsCoveragePerEvidence[annotationsCoveragePerEvidence$evidence == 'curation_effort',]
byCuration_effortDF$propor_TargetsAnnotationsInEvidence

# byN_references NO APORTA NADA
#byN_references <- annotationsCoveragePerEvidence[annotationsCoveragePerEvidence$evidence == 'n_references',]

ggplot(byCuration_effortDF, 
       aes(x=minEvidenceScore, y=propor_TargetsAnnotationsInEvidence,
           color=annotation, group=annotation)) +
  scale_color_manual(values = color_vals2, breaks = names(color_vals2)) +
  geom_line(size=1, linetype='solid') +
  theme(panel.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.2,"cm"),
        legend.key = element_rect(fill = "white"),
        plot.margin = margin(l = 0, unit = "cm"))+
  ylab("Annotation coverage (%)")+
  xlab("Evidence Level") +
  labs(fill = "Annotation")



a <- ggplot(regulons_coverage, aes(x = coverage, y = database))+
  geom_col(aes(fill = interaction), position = "dodge2", color = "black", linewidth = 0.1)+
  scale_fill_manual(values = color_vals, breaks = names(color_vals))+
  scale_x_continuous(position = "top")+
  theme(panel.background = element_blank(),
        axis.title = element_text(size = 4),
        axis.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.2,"cm"),
        legend.key = element_rect(fill = "white"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(l = -0.2, unit = "cm"))+
  xlab("Annotation coverage (%)")+
  labs(fill = "Annotation")



View(annotationsCoveragePerEvidence)

annotationsCoveragePerEvidence$
  
  
  regulons_coverage <- rbindlist(lapply(annotations, function(ann){
    ann_table <- ann_tables %>% filter(annotation == ann)
    
    ## Number Gene-Term Pairs
    n <- nrow(unique(ann_table))
    
    confidence_level <- rbindlist(lapply(c("A","B","C","D","E"), function(conf){
      regulons_filt <- regulons %>% filter(confidence <= conf) %>% unique()
      target_genes <- regulons_filt %>% select(target) %>% unique() %>% pull(target)
      tfs <-  regulons_filt %>% select(tf) %>% unique() %>% pull(tf)
      ann_table_target <- ann_table %>% filter(symbol %in% target_genes) %>% unique()
      ann_table_tfs <- ann_table %>% filter(symbol %in% tfs) %>% unique()
      
      coverage_target <- round((nrow(ann_table_target)/n)*100,2)
      coverage_tf <- round((nrow(ann_table_tfs)/n)*100,2)
      stats <- data.frame(coverage = c(coverage_target, coverage_tf), type = c("Target Genes","TFs"), database = conf,annotation = ann) %>%
        mutate(interaction = paste0(annotation,"\nwith ",type))
    }))
    return(confidence_level)
  })) %>% mutate(interaction = factor(interaction, levels = unique(interaction)))

color_vals <- c()
for (annotation in annotations){
  base_color <- colors_blind[names(colors_blind) %in% annotation]
  dark_color <- darken(base_color, 0.3)
  light_color <- lighten(base_color,0.3)
  color_vals_insert <- c(dark_color, light_color)
  names(color_vals_insert) <- c(paste0(annotation,"\nwith Target Genes"),
                                paste0(annotation,"\nwith TFs"))
  color_vals <- c(color_vals, color_vals_insert)
}


evidences <- data.frame(database = c("A","A","A","A","A","B","B","B","B","C","C","C","D","D","E","E"),
                        evidence = c("> 2 curated\nliterature", "1 curated\nliterature", "ChIP-Seq","TFBMs Prediction","GTEx inference",
                                     "1 curated\nliterature","ChIP-Seq","TFBMs Prediction","GTEx inference",
                                     "1 curated\nliterature","ChIP-Seq","TFBMs Prediction",
                                     "1 curated\nliterature","ChIP-Seq",
                                     "TFBMs Prediction","GTEx inference")) %>%
  mutate(evidence = factor(evidence, levels = unique(evidence)),
         database = factor(database, levels = rev(unique(database))))

color_evidences <- c("#c7d8e0","#c7e1d4","#f1d0ad","#f1e4af","#f1e4af")
names(color_evidences) <- levels(evidences$evidence)

regulons_coverage$interaction <- factor(regulons_coverage$interaction, levels = rev(levels(regulons_coverage$interaction)))
regulons_coverage$database <- factor(regulons_coverage$database, levels = levels(evidences$database))

a <- ggplot(regulons_coverage, aes(x = coverage, y = database))+
  geom_col(aes(fill = interaction), position = "dodge2", color = "black", linewidth = 0.1)+
  scale_fill_manual(values = color_vals, breaks = names(color_vals))+
  scale_x_continuous(position = "top")+
  theme(panel.background = element_blank(),
        axis.title = element_text(size = 4),
        axis.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.2,"cm"),
        legend.key = element_rect(fill = "white"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(l = -0.2, unit = "cm"))+
  xlab("Annotation coverage (%)")+
  labs(fill = "Annotation")

b <- ggplot(evidences, aes(x = evidence, y = database))+
  geom_point(aes(fill = evidence), fill = "gray80", size = 1, color = "black", shape = 21, stroke = 0.1)+
  scale_x_discrete(position = "top")+
  theme(panel.background = element_blank(),
        axis.title = element_text(size = 4),
        axis.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.2,"cm"),
        legend.key = element_rect(fill = "white"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45),
        axis.ticks = element_blank())+
  xlab("Evidences")+
  ylab("Regulon Confidence Levels")+
  labs(fill = "Evidences")

gg_2 <- plot_grid(b,a, ncol = 2, align = "h",
                  rel_widths = c(.5,1))

ggsave("Figures/Figure 2 v2.tiff", plot = gg_2,units = "cm",height = 6, width = 8,dpi = 500,
       compression = "lzw",bg = "white")

dorothea_info <- vroom("data/dorothea.tsv") %>%
  as.data.frame() %>%
  filter(org == 9606, confidence <= "C") %>%
  select(tf,target,confidence)
target_genes <- unique(dorothea_info$target)

count_genes <- do.call("rbind",lapply(annotations, function(annotation){
  data_content <- ann_tables %>% filter(annotation == !!annotation, symbol %in% target_genes)
  tfs_by_target <- dorothea_info %>% filter(target %in% unique(data_content$symbol)) %>%
    group_by(target) %>%
    summarise(tfs = n()) %>%
    mutate(annotation = annotation)
  return(tfs_by_target)
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
  count_genes_top_table <- count_genes_top %>% filter(annotation == !!annotation)
  ann_table <- ann_tables %>% filter(annotation == !!annotation, symbol %in% count_genes_top_table$target)
  ann_table_show <- as.data.frame(table(ann_table$annotation_id)) %>%
    arrange(desc(Freq)) %>%
    do(head(.,10))
  ann_table_show$annotation_id <- as.character(ann_table_show$Var1)
  ann_table_show$annotation <- annotation
  ann_table_show <- merge(ann_table_show,ann_info,by = "annotation_id") %>%
    select(annotation_id,annotation,term,Freq)
  
  # to_file <- ann_table_show[order(ann_table_show$Freq,decreasing = T),c("annotation_id","term","Freq")] %>%
  #   as.data.frame()
  # colnames(to_file) <- c("Annotation Id","Term","Number of top genes implicated in gene set")
  # 
  # options(java.parameters = "-Xmx8000m")
  # message(glue("Writing table for {annotation}"))
  # if (annotation == "KEGG"){
  #   invisible(file.remove("Tables/Supplementary Table 4.xlsx"))
  #   write.xlsx(to_file, file = "Tables/Supplementary Table 4.xlsx",
  #              append = FALSE,
  #              sheetName = annotation,
  #              row.names = F)
  # } else{
  #   write.xlsx(to_file, file = "Tables/Supplementary Table 4.xlsx",
  #              append = T,
  #              sheetName = annotation,
  #              row.names = F)
  # }
  return(ann_table_show)
}))

ann_table <- ann_table[order(ann_table$Freq),]
ann_table$annotation_id <- factor(ann_table$annotation_id,levels = unique(ann_table$annotation_id))
ann_table$annotation <- factor(ann_table$annotation, levels = annotations)

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


########## Null simulations results ##########

list_files <- list.files("random_list_analysis",pattern = "*tsv",full.names = T)

fullContent <- rbindlist(lapply(list_files, function(fileName){
  fileNameCheck <- gsub("GO_BP","GO BP",fileName)
  basenameFile <- basename(fileNameCheck)
  annotation <- unlist(strsplit(basenameFile,"_"))[1]
  size = as.numeric(gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[2]))
  content <- vroom(fileName) %>% mutate(annotation = annotation, size = size, interaction = paste0(typeOfAnalysis,"-",annotation))
})) %>% mutate(typeOfAnalysis = factor(typeOfAnalysis, levels = c("Target_NonCentral","Target_Hypergeom")),
               size = factor(size, levels = unique(size[order(size)])),
               annotation = factor(annotation, levels = unique(annotation)))

fullContent <- fullContent %>% filter(typeOfAnalysis %in% c("Target_Hypergeom","Target_NonCentral"))

shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}

annotation = "GO BP"
fullPlots <- lapply(annotations, function(annotation){
  base_color <- colors_blind[names(colors_blind) %in% annotation]
  dark_color <- darken(base_color, 0.3)
  light_color <- lighten(base_color,0.3)
  color_vals_insert <- c(dark_color, light_color)
  names(color_vals_insert) <- c("Fisher's Exact Test","Wallenius Test")
  
  fullContentAnn <- fullContent %>% filter(annotation == !!annotation, size %in% c(2,5,10,20,30)) %>%
    arrange(size) %>%
    mutate(typeOfAnalysis = ifelse(typeOfAnalysis == "Target_Hypergeom","Fisher's Exact Test","Wallenius Test"),
           typeOfAnalysis = factor(typeOfAnalysis, levels = c("Fisher's Exact Test","Wallenius Test")),
           sizeName = factor(paste0("TFs Size ",size), levels = unique(paste0("TFs Size ",size))),
           percentage = p_times * 100)
  
  hist_info <- rbindlist(lapply(levels(fullContentAnn$sizeName), function(size){
    hist_info <- fullContentAnn %>% filter(sizeName == !!size)
    hist_info <- rbindlist(lapply(levels(hist_info$typeOfAnalysis), function(typeOfAnalysis){
      hist_info <- hist_info %>% filter(typeOfAnalysis == !!typeOfAnalysis)
      p <- ggplot(hist_info, aes(x = percentage))+
        geom_histogram(breaks = seq(0,100,by = 2))
      hist_info <- ggplot_build(p)$data[[1]]
      data.frame(counts = hist_info$count,percentage = hist_info$x, xmin = hist_info$xmin, xmax = hist_info$xmax,
                 typeOfAnalysis = typeOfAnalysis, sizeName = size) %>% filter(counts != 0)
    }))
  })) %>% mutate(typeOfAnalysis = factor(typeOfAnalysis, levels = c("Fisher's Exact Test","Wallenius Test")),
                 sizeName = factor(sizeName, levels = levels(fullContentAnn$sizeName)))
  
  p <- ggplot(hist_info)+
    facet_wrap(~sizeName, ncol = 2, scales = "free")+
    geom_col(aes(x = percentage, y = counts, fill = typeOfAnalysis),linewidth = 0.2, position = "identity", color = "black", alpha = 0.5, width = 2)+
    #geom_vline(xintercept = 0.05, linewidth = 0.2, color = "red", linetype = "dotted")+
    scale_fill_manual(values = color_vals_insert)+
    theme(panel.background = element_blank(),
          strip.text = element_text(size = 4),
          strip.background = element_rect(fill = "white"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 4),
          title = element_text(face = "bold"),
          axis.title = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 4),
          legend.key.size = unit(0.2,"cm"))+
    ylab("Frequency")+
    xlab("Percentage of times significantly")+
    scale_y_log10()+
    labs(fill = annotation)+
    scale_x_continuous(limits = c(0,100))
    
  
  p <- shift_legend2(p)
  dev.off()
  
  content_top <- fullContentAnn %>% filter(annotation_id %in% ann_table$annotation_id)
  
  p2 <- ggplot(content_top, aes(x = size, y = ranking))+
    geom_jitter(aes(color = typeOfAnalysis), size = 0.3, width = 0.2)+
    scale_color_manual(values = color_vals_insert)+
    theme(panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 4),
          title = element_text(face = "bold"),
          strip.text = element_text(size = 4,face = "bold"),
          axis.title = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          legend.margin=margin(0,8,0,0),
          legend.box.margin=margin(-10,0,-10,-10),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          legend.position = "none",
          plot.title = element_text(size = 4, hjust = 0.5))+
    ggtitle(annotation)+
    xlab("TFs Size")+
    ylab("Median Rank of Top10 Terms")
  
  return(list(p,p2))
})

gg1 <- plot_grid(plotlist = sapply(fullPlots,"[",1))
gg2 <- plot_grid(plotlist = sapply(fullPlots,"[",2), nrow = 1)

gg4_v2 <- ggdraw()+
  draw_plot(gg1,x = 0, y = 0.2, height = 0.8, width = 1)+
  draw_plot(gg2, x = 0, y = 0, height = 0.2, width = 1) +
  draw_plot_label(c("A","B"), x = c(0, 0), y = c(1, 0.225), fontface = "plain", family = "serif", size = 10)

# gg4_v2 <- plot_grid(gg4_v2, gg_legend, rel_widths = c(0.9,0.2))
ggsave("Figures/Figure 4 v2.tiff", plot = gg4_v2,
       units = "cm",height = 18, width = 20,dpi = 500,
       compression = "lzw", bg = "white")


list_files <- list.files("random_list_analysis",pattern = "*tsv",full.names = T)

fullContent <- rbindlist(lapply(list_files, function(fileName){
  fileNameCheck <- gsub("GO_BP","GO BP",fileName)
  basenameFile <- basename(fileNameCheck)
  annotation <- unlist(strsplit(basenameFile,"_"))[1]
  size = as.numeric(gsub(".tsv","",unlist(strsplit(basenameFile,"_"))[2]))
  content <- vroom(fileName) %>% mutate(annotation = annotation, size = size, interaction = paste0(typeOfAnalysis,"-",annotation))
})) %>% mutate(size = factor(size, levels = unique(size[order(size)])),
               annotation = factor(annotation, levels = unique(annotation)))

fullContent <- fullContent %>% filter(typeOfAnalysis %in% c("TFs_Hypergeom"))

annotation = "GO BP"
fullPlots <- lapply(annotations, function(annotation){
  base_color <- colors_blind[names(colors_blind) %in% annotation]
  names(base_color) <- NULL
  
  fullContentAnn <- fullContent %>% filter(annotation == !!annotation, size %in% c(2,5,10,20,30)) %>%
    arrange(size) %>%
    mutate(typeOfAnalysis = "Fisher's Exact Test",
           sizeName = factor(paste0("TFs Size ",size), levels = unique(paste0("TFs Size ",size))),
           percentage = p_times * 100)
  
  hist_info <- rbindlist(lapply(levels(fullContentAnn$sizeName), function(size){
    hist_info <- fullContentAnn %>% filter(sizeName == !!size)
    hist_info <- rbindlist(lapply(unique(hist_info$typeOfAnalysis), function(typeOfAnalysis){
      hist_info <- hist_info %>% filter(typeOfAnalysis == !!typeOfAnalysis)
      p <- ggplot(hist_info, aes(x = percentage))+
        geom_histogram(breaks = seq(0,100,by = 2))
      hist_info <- ggplot_build(p)$data[[1]]
      data.frame(counts = hist_info$count,percentage = hist_info$x, xmin = hist_info$xmin, xmax = hist_info$xmax,
                 typeOfAnalysis = typeOfAnalysis, sizeName = size) %>% filter(counts != 0)
    }))
  })) %>% mutate(sizeName = factor(sizeName, levels = levels(fullContentAnn$sizeName)))
  
  p <- ggplot(hist_info)+
    facet_wrap(~sizeName, ncol = 2, scales = "free")+
    geom_col(aes(x = percentage, y = counts, fill = typeOfAnalysis),linewidth = 0.2, position = "identity", color = "black", alpha = 0.5, width = 2)+
    #geom_vline(xintercept = 0.05, linewidth = 0.2, color = "red", linetype = "dotted")+
    scale_fill_manual(values = base_color)+
    theme(panel.background = element_blank(),
          strip.text = element_text(size = 4),
          strip.background = element_rect(fill = "white"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 4),
          title = element_text(face = "bold"),
          axis.title = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 4),
          legend.position = "none",
          legend.key.size = unit(0.2,"cm"),
          plot.title = element_text(size = 4, hjust = 0.5))+
    ggtitle(annotation)+
    ylab("Frequency")+
    xlab("Percentage of times significantly")+
    scale_y_log10()+
    labs(fill = annotation)+
    scale_x_continuous(limits = c(0,100))
  
  content_top <- fullContentAnn %>% arrange(desc(p_times)) %>% pull(annotation_id) %>% unique() %>% head(6)
  content_top <- fullContentAnn %>% filter(annotation_id %in% content_top) %>% inner_join(ann_info) %>% arrange(desc(p_times)) %>%
    mutate(term = factor(term, levels = unique(term)),
           annotation_id = factor(annotation_id, levels = unique(annotation_id)))
  
  p2 <- ggplot(content_top, aes(x = size, y = percentage))+
    facet_wrap(~annotation_id, scales = "free")+
    geom_col(fill = base_color)+
    theme(panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 4),
          title = element_text(face = "bold"),
          strip.text = element_text(size = 3,face = "bold"),
          axis.title = element_text(size = 5),
          axis.text.y = element_text(size = 4),
          axis.line.x = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth=0.2, linetype='solid'),
          legend.margin=margin(0,8,0,0),
          legend.box.margin=margin(-10,0,-10,-10),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          legend.position = "none",
          plot.title = element_text(size = 4, hjust = 0.5))+
    ggtitle(annotation)+
    xlab("TFs Size")+
    ylab("Percentage of times significantly")+
    scale_y_continuous(limits = c(0,100))

  return(list(p,p2))
})

gg1 <- plot_grid(plotlist = sapply(fullPlots,"[",1))
gg2 <- plot_grid(plotlist = sapply(fullPlots,"[",2), nrow = 1)

gg4_v2 <- ggdraw()+
  draw_plot(gg1,x = 0, y = 0.3, height = 0.7, width = 1)+
  draw_plot(gg2, x = 0, y = 0, height = 0.3, width = 1) +
  draw_plot_label(c("A","B"), x = c(0, 0), y = c(1, 0.325), fontface = "plain", family = "serif", size = 10)

# gg4_v2 <- plot_grid(gg4_v2, gg_legend, rel_widths = c(0.9,0.2))
ggsave("Figures/Supplementary Figure 1 v2.tiff", plot = gg4_v2,
       units = "cm",height = 18, width = 20,dpi = 500,
       compression = "lzw", bg = "white")

