library(vroom)
library(data.table)
library(tidyverse)
library(parallel)
library(glue)
library(lubridate)

# PID 3486739 <- 

resultsFiles <- Sys.glob('random_lists_analysis/TFsLists/*/*/[0-9]*Merged*.txt', dirmark = FALSE)
annotNsize <- do.call(rbind, strsplit(resultsFiles,split = '/'))[,3:4]
annotNsize <- as.data.frame(annotNsize[!duplicated(annotNsize),])
colnames(annotNsize) <- c("db","size")

not_done <- nrow(annotNsize)

masterSummaryDF <- lapply(1:nrow(annotNsize), function(idx){
  st <- Sys.time()
  message(glue("Summary Random Simulations Results [{idx}/{nrow(annotNsize)}]"))
  db <- annotNsize$db[idx]
  size <- annotNsize$size[idx]
  resultsFiles <- Sys.glob(file.path('random_lists_analysis/TFsLists',db,size,'[0-9]*Merged*.txt'), dirmark = FALSE)
  size <- as.numeric(gsub("size","",size))
  message(glue("\tSummary of {db} with size {size}"))
  chunks <- split(resultsFiles, ceiling(seq_along(resultsFiles)/15))
  statisticsOutFile <- file.path(dirname(resultsFiles[1]),'_statistics.tsv')
  not_done <<- not_done - 1
  if (!file.exists(statisticsOutFile)){
    message("\t\tCalculating Ranks")
    masterDF <- rbindlist(mclapply(chunks, function(x_sub){
      masterDF <- rbindlist(lapply(x_sub, function(resultsFile){
        results <- vroom(resultsFile)
        results <- results %>%
          mutate(rankTFs = rank(pvals_TFsFisher),
                 rankTargetsF = rank(pvals_targetsFisher),
                 rankTargetsW = rank(pvals_targetsWallenius))
        results <- results %>% select(term, annotation, k, pvals_TFsFisher, rankTFs, pvals_targetsFisher, rankTargetsF, pvals_targetsWallenius, rankTargetsW)
      }))
      return(masterDF)
    }, mc.cores = 15))
    message("\t\tCalculating Stats")
    statsMasterDF <- masterDF %>% group_by(term) %>% # group by term
      summarise(pTFs = mean(pvals_TFsFisher < 0.05), # proportion of p < 0.05 for TFs by Fisher
                rTFs = median(rankTFs), # median ranking according to TFs by Fisher
                Q1rTFs = quantile(rankTFs, 0.25), # Q1 ranking according to TFs by Fisher
                Q3rTFs = quantile(rankTFs, 0.75), # Q3 ranking according to TFs by Fisher
                pTargetsF = mean(pvals_targetsFisher < 0.05), # proportion of p < 0.05 for Targets by Fisher
                rTargetsF = median(rankTargetsF), # median ranking according to Targets by Fisher
                Q1rTargetsF = quantile(rankTargetsF, 0.25), # Q1 ranking according to Targets by Fisher
                Q3rTargetsF = quantile(rankTargetsF, 0.75), # Q3 ranking according to Targets by Fisher
                pTargetsW = mean(pvals_targetsWallenius < 0.05), # proportion of p < 0.05 for Targets by Wallenius
                rTargetsW = median(rankTargetsW), # median ranking according to Targets by Wallenius
                Q1rTargetsW = quantile(rankTargetsW, 0.25), # Q1 ranking according to Targets by Wallenius
                Q3rTargetsW = quantile(rankTargetsW, 0.75)) %>% # Q3 ranking according to Targets by Wallenius
      mutate(annotation = db, size = size)
    statsMasterDF <- statsMasterDF %>% mutate(rTFs = ifelse(is.na(pTFs),NA,rTFs),
                                              Q1rTFs = ifelse(is.na(pTFs),NA,Q1rTFs),
                                              Q3rTFs = ifelse(is.na(pTFs),NA,Q3rTFs),
                                              rTargetsF = ifelse(is.na(pTargetsF),NA,rTargetsF),
                                              Q1rTargetsF = ifelse(is.na(pTargetsF),NA,Q1rTargetsF),
                                              Q3rTargetsF = ifelse(is.na(pTargetsF),NA,Q3rTargetsF),
                                              rTargetsW = ifelse(is.na(pTargetsW),NA,rTargetsW),
                                              Q1rTargetsW = ifelse(is.na(pTargetsW),NA,Q1rTargetsW),
                                              Q3rTargetsW = ifelse(is.na(pTargetsW),NA,Q3rTargetsW))
    
    message(glue("\t\tWriting file {statisticsOutFile}"))
    write.table(statsMasterDF, statisticsOutFile, sep='\t',quote = F,row.names = F,col.names = T)
    et <- Sys.time()
    seconds <- round(as.numeric(difftime(time1 = et, time2 = st, units = "secs")), 3)
    estimated_time <- seconds * not_done
    message(glue("\tSummary of {db} with size {size} DONE"))
    estimated_time <- seconds_to_period(estimated_time)
    message(glue("{estimated_time} to finish all analysis"))
  }
  return(statisticsOutFile)
})

masterSummaryDF <- rbindlist(lapply(masterSummaryDF, function(x){vroom(x)}))

masterSummaryOutFile <- 'random_lists_analysis/masterSimsResults.tsv'
write.table(masterSummaryDF, masterSummaryOutFile, sep='\t',quote = F,row.names = F,col.names = T)



