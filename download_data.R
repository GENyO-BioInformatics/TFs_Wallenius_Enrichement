library(optparse)
library(glue)


option_list <- list(make_option(c("-o", "--organism"),
    type = "character", default = "9606",
    help = "organism code", metavar = "character"
))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
organism <- opt$organism


##############################
## CollecTRI
##############################


## Fetch Data Clean and Save TF-Target Data
# dorotheaTFsGRN <- OmnipathR::dorothea(organism=9606, genesymbols=TRUE, loops=TRUE)
collectriTFsGRN <- OmnipathR::collectri(organism = organism, genesymbols = TRUE, loops = TRUE)
collectriTFsGRN <- as.data.frame(collectriTFsGRN)

# collectriTFsGRN$source_genesymbol ## TF Gene
# collectriTFsGRN$target_genesymbol ## Target Gene
collectriTFsGRN_clean <- data.frame(tf = collectriTFsGRN$source_genesymbol, target = collectriTFsGRN$target_genesymbol, confidence = "A", org = "9606")

write.table(collectriTFsGRN, file = glue("data/collectri_raw_{organism}.tsv"), sep = "\t", row.names = F, col.names = T)
write.table(collectriTFsGRN_clean, file = glue("data/collectri_{organism}.tsv"), sep = "\t", row.names = F, col.names = T)


##############################
## Annotation Databases
##############################

annotations <- c("KEGG","GO_BP","WikiPathways","Reactome")
endpoint <- "https://genecodis.genyo.es/gc4/database?"
organism <- "9606"
for (annotation in annotations){
  query <- paste0(endpoint,"annotation=",annotation,"&nomenclature=symbol&organism=",organism)
  print(query)
  system(paste0("wget '",query,"' -O ",annotation,".tsv"))
}
