library(vroom)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(colorspace)
library(glue)
library(data.table)
library(parallel)
colors_blind <- c("KEGG" = "#979A61",
                  "GO BP" = "#3d91e0",
                  "Reactome" = "#f1b620",
                  "WikiPathways" = "#D974A0")

# cargar collectri
collectri <- read.delim("../collectri.tsv")
annotations <- names(colors_blind)

set.seed(123456789)
ann_prob <- lapply(annotations, function(ann){
  print(ann)
  anname <- gsub(" ","_",ann)
  # cargar base de datos
  anntable <- read.delim(glue("../{anname}.tsv")) %>% filter(organism == 9606)
  # Filtrar con genes anotados en base de datos y collectri
  collectri_ann <- collectri %>% filter(target %in% anntable$symbol)
  anntable <- anntable %>% filter(symbol %in% collectri_ann$target)
  ids <- unique(anntable$annotation_id)
  # número de tfs
  num_tfs <- length(unique(collectri_ann$tf))
  # dataframe con tf y lista de target genes
  targets_per_tf <- collectri_ann %>%
    filter(target %in% anntable$symbol) %>%
    group_by(tf) %>%
    summarize(targets = list(target))
  monte_carlo_iterations <- 10
  
  # Pre-compute unique TFs and their target genes
  unique_tfs <- unique(collectri_ann$tf)
  tf_targets <- lapply(unique_tfs, function(u_tf) {
    collectri_ann %>% filter(tf == u_tf) %>% pull(target)
  })
  # Create a named list for quick lookup
  names(tf_targets) <- unique_tfs
  
  # Function to check if id is in the annotation table for a given TF
  check_id_in_annotations <- function(id, tf_targets, anntable) {
    tf_random <- sample(names(tf_targets), size = 1)
    targets <- tf_targets[[tf_random]]
    anntable_targets <- anntable %>% filter(symbol %in% targets)
    return(id %in% anntable_targets$annotation_id)
  }
  
  # Main computation
  sim_probs <- mclapply(ids, function(id) {
    results <- lapply(1:monte_carlo_iterations, function(i) {
      check_id_in_annotations(id, tf_targets, anntable)
    })
    mean(unlist(results))
  }, mc.cores = 10)
  
  probs <- mclapply(ids, function(id){
    # que genes están en la anotacion?
    genes_with_annotation <- anntable %>%
      filter(annotation_id == id) %>%
      pull(symbol) %>%
      unique()
    # cuantos de estos genes son targets de TFs?
    # dataframe con TFs y count_with_annotation
    targets_with_annotation <- targets_per_tf %>%
      rowwise() %>%
      mutate(count_with_annotation = sum(targets %in% genes_with_annotation))
    # cuantos TFs tienen al menos 1 gen en la anotación?
    num_tfs_with_annotation <- sum(targets_with_annotation$count_with_annotation > 0)
    probability <- num_tfs_with_annotation / num_tfs
    return(probability)
  }, mc.cores = 10)
  probs <- data.frame(sims = unlist(sim_probs), prob = unlist(probs), annotation_id = ids)
  probs <- probs %>% rename(probability = prob, simulations = sims) %>% mutate(ann = ann) %>%
    dplyr::select(ann, annotation_id, probability, simulations)
  return(probs)
})

ann_plots <- lapply(ann_prob, function(probs){
  max_prob <- probs %>% top_n(n = 5, wt = probability) %>% arrange(desc(probability))
  brks <- seq(min(probs$probability), max(probs$probability), length.out = 30)
  H <- hist(probs$probability, plot = FALSE, breaks = brks)
  # segment
  y <- max(H$counts) / 10
  x <- 0.8
  y_0 <- max(H$counts) / 6.5
  y_6 <- max(H$counts) / 1.5
  ys <- seq(y_6, y_0, length.out = 6)
  
  p <- ggplot()+
    geom_histogram(data = probs, aes(x = probability, fill = ann), breaks = brks, color = "gray20", linewidth = 0.2)+
    theme_classic(base_size = 4,base_line_size = 0.2)+
    theme(legend.position = "none",
          strip.background = element_blank())+
    ylab("Frequency")+
    xlab("Probabilities of Annotation IDs")+
    annotate("text", x = x, y = ys[1], label = "Highest probabilities", size = 3/.pt)+
    scale_fill_manual(values = colors_blind)+
    facet_wrap(~ann)+
    # geom_segment(data = data.frame(xmin = xmin, 
    #                                xmax = max(max_prob$probability), 
    #                                ymin = y, 
    #                                ymax = y),
    #              aes(x = xmin, xend = xmax, y=ymin, yend=ymax),
    #              linewidth = 0.2)+
    # geom_segment(data = data.frame(xmin = xmin, 
    #                                xmax = xmin, 
    #                                ymin = y, 
    #                                ymax = y*0.6),
    #              aes(x = xmin, xend = xmax, y=ymin, yend=ymax),
    #              linewidth = 0.2)+
    # geom_segment(data = data.frame(xmin = max(max_prob$probability), 
    #                                xmax = max(max_prob$probability), 
    #                                ymin = y, 
    #                                ymax = y*0.6),
    #              aes(x = xmin, xend = xmax, y=ymin, yend=ymax),
    #              linewidth = 0.2)+
    scale_x_continuous(limits = c(0, 1))
  for (i in 1:nrow(max_prob)){
    p <- p + annotate("text", x = x, y = ys[i+1], label = paste0(max_prob$annotation_id[i],": ",round(max_prob$probability[i],3)),  size = 3/.pt)
  }
  p
  return(p)
})

ann_prob_tbl <- rbindlist(ann_prob)
write.csv(ann_prob_tbl, file = "probabilities.csv", row.names = F, quote = F)


library(cowplot)
p <- cowplot::plot_grid(plotlist = ann_plots)

ggsave("figure_probabilities.tiff", plot = p,units = "cm",height = 6, width = 8,dpi = 500,
       compression = "lzw",bg = "white")

ggsave("figure_probabilities.png", plot = p,units = "cm",height = 6, width = 8,dpi = 500,
       bg = "white")
