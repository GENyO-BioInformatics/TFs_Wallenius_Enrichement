# Load required libraries
library(parallel)
library(data.table)
library(tidyverse)
library(glue)
library(reshape)
library(pbapply)

# Function to analyze a single TF annotation file
analyze_file <- function(f, ann) {
  require(glue) # Load the glue package for string interpolation
  require(tidyverse) # Load the tidyverse package for data manipulation
  require(reshape) # Load reshape package
  
  # Extract the name prefix from the file path
  n <- basename(f)
  n <- unlist(strsplit(n, "_"))[1]
  
  # Generate the output file path for the TF enrichment results
  outfile_tfs <- gsub("TFsAnnoted", "TFsAnnotedRes", f)
  
  # If the TF enrichment results do not exist, run the Python enrichment script
  if (!file.exists(outfile_tfs)) {
    system(
      glue::glue(
        "python TFsEnrichment.py -ann {ann} -tfs {f} -out {outfile_tfs} -m tfs"
      )
    )
  }
  
  # Update the file path to the annotated targets file
  f <- gsub("TFsAnnoted", "TFsTargetsAnnoted", f)
  outfile_targetsF <- gsub("TFsTargetsAnnoted", "TFsTargetsFAnnotedRes", f)
  
  # Targets Fisher
  if (!file.exists(outfile_targetsF)) {
    system(
      glue::glue(
        "python TFsEnrichment.py -ann {ann} -tfs {f} -out {outfile_targetsF} -m targetsF"
      )
    )
  }
  
  outfile_targetsW <- gsub("TFsTargetsAnnoted", "TFsTargetsWAnnotedRes", f)
  # Targets Wallenius
  if (!file.exists(outfile_targetsW)) {
    system(
      glue::glue(
        "python TFsEnrichment.py -ann {ann} -tfs {f} -out {outfile_targetsW} -m targetsW"
      )
    )
  }
  
  k <- dirname(f)
  k <- basename(k)
  k <- as.numeric(gsub("size","",k))
  
  # Read and select relevant columns from the TF enrichment results
  outfile_tfs <- read.delim(outfile_tfs) %>% select(term, pvals_TFsFisher)
  
  # Read and select relevant columns from Target Fisher
  outfile_targetsF <- read.delim(outfile_targetsF) %>% select(term, pvals_targetsFisher)
  
  # Read and select relevant columns from the target enrichment results
  outfile_targetsW <- read.delim(outfile_targetsW) %>% select(term, pvals_targetsWallenius)
  
  # merge three datasets
  df_list <- list(outfile_tfs, outfile_targetsF, outfile_targetsW)
  merged <- merge_recurse(df_list) %>% mutate(annotation = ann, k = k)
  
  # write table
  outfile <- gsub("TFsTargetsAnnoted","MergedRes",f)
  if (!file.exists(outfile)){
    write.table(merged, file = outfile, quote = F, row.names = F, sep = "\t")
  }
}


# List all annotation folders in the specified directory
annotations = list.files(path = "random_lists_analysis/TFsLists/")

# Iterate over each annotation file
res <- lapply(annotations, function(ann) {
  message(glue("Analyse {ann} ...")) # Print progress message
  
  # List all size-specific directories within the annotation directory
  sizes <- list.files(file.path("random_lists_analysis/TFsLists", ann))
  
  # Iterate over each size-specific directory
  results <- lapply(sizes, function(s) {
    message(glue("\tSize {s}")) # Print progress message
    
    # List all TF annotation files for the current size
    x <- list.files(file.path("random_lists_analysis/TFsLists", ann, s),
                    pattern = "TFsAnnoted.txt",
                    full.names = T)
    
    # check if final files exist
    end_files <- gsub("TFsAnnoted","MergedRes",x)
    if (!all(file.exists(end_files))){ # if not: we need to launch code to each file
      # divide 1000 files in 10 chunks of 100 files each one
      chunks <- split(x, rep_len(1:10, length(x)))
      # each chunk will be processed in a lapply
      # this will show a progress bar of each 100 files
      res <- pblapply(chunks, function(y){
        # divide 100 files in 15 chunks
        sub_chunks <- split(y, rep_len(1:15, length(y)))
        # launch each sub chunk with mclapply
        res <- mclapply(sub_chunks, function(z){
          res <- lapply(z, analyze_file, ann = ann)
        }, mc.cores = 15)
      })
    }
  })
})
