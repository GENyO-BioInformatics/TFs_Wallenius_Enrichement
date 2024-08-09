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
  content_top <- fullContentAnn %>% filter(annotation_id %in% content_top) %>% arrange(desc(p_times)) %>%
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
ggsave("MegaFigure.tiff", plot = gg4_v2,
       units = "cm",height = 18, width = 20,dpi = 500,
       compression = "lzw", bg = "white")

