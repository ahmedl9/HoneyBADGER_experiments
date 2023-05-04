#In case of corruption of library can delete before adding back
#rm(list = ls(all = TRUE))
#remove.packages("HoneyBADGER")

#Please change wd to wherever project is stored. 
setwd("C:/Users/Ahmed Lone/Desktop/COMS/final_project/rstudio")
devtools::install("HoneyBADGER")
library(HoneyBADGER)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("dbscan", quietly = TRUE)) {
  install.packages("dbscan")
}

library(dplyr)
library(reshape2)


# Function to simulate deletions
simulate_deletions <- function(r, deletion_size, clonality) {
  r_copy <- r
  num_cells <- ncol(r_copy)
  num_deletions <- 5
  num_rows <- nrow(r_copy)
  
  # Calculate number of cells affected by deletions
  num_affected_cells <- ceiling(num_cells * clonality)
  
  # Generate deletion locations
  deletion_starts <- sample(1:(num_rows - deletion_size), num_deletions, replace = TRUE)
  
  # Modify r_copy matrix
  for (i in 1:num_deletions) {
    affected_cells <- sample(1:num_cells, num_affected_cells)
    r_copy[(deletion_starts[i]:(deletion_starts[i] + deletion_size - 1)), affected_cells] <- 0
  }
  
  return(r_copy)
}


# Function to compute sensitivity and precision of the honeybadger model
calculate_sensitivity_precision <- function(r, r_simulated, binary_matrix) {
  true_positives <- sum((binary_matrix == 1) & (r_simulated == r))
  false_positives <- sum((binary_matrix == 1) & (r_simulated != r))
  false_negatives <- sum((binary_matrix == 0) & (r_simulated != r))
  
  sensitivity <- true_positives / (true_positives + false_negatives)
  precision <- true_positives / (true_positives + false_positives)
  
  return(list(sensitivity = sensitivity, precision = precision))
}

process_results <- function(r, r_simulated, results) {
  # Get real deletions by assessing simulated matrix. Find cell/chromosome combinations of real matrix 
  r_long <- melt(r, variable.name = "sample", value.name = "r_value")
  r_mod_long <- melt(r_simulated, variable.name = "sample", value.name = "r_mod_value")
  combined_data <- merge(r_long, r_mod_long, by = c("Var1", "Var2"))
  real_deletions <- combined_data #subset(combined_data, combined_data$r_value != combined_data$r_mod_value)
  
  real_deletions$chr <- sapply(strsplit(as.character(real_deletions$Var1), split = ":"), "[", 1)
  real_deletions$deletion <- real_deletions$r_value != real_deletions$r_mod_value
  real_deletions <- real_deletions[, !colnames(real_deletions) %in% c("Var1")]
  
  real_deletions_by_chr <- group_by(real_deletions, chr, Var2)
  #real_deletions_by_chr <- slice_min(real_deletions_by_chr, order_by = row_number(), n = 1)
  real_deletions_by_chr <- summarize_all(real_deletions_by_chr, funs(sum(., na.rm = TRUE)))
  real_deletions_by_chr <- ungroup(real_deletions_by_chr)
  
  #add chr number to results
  results_long <- results[, !colnames(results) %in% c("start", "end", "width", "strand", "avg.del.loh.allele")]
  results_long <- melt(results_long, id.vars = "seqnames", variable.name = "Var2", value.name = "detected_value")
  results_long$chr <- sapply(strsplit(as.character(results_long$seqnames), split = "chr"), "[", 2)
  
  final_results <- merge(real_deletions_by_chr, results_long, by = c("chr", "Var2"), all.x = TRUE)
  final_results[is.na(final_results)] <- 0
  
  final_results$tp = (final_results$r_value != final_results$r_mod_value) & (final_results$detected_value >= 0.5)
  final_results$fp = (final_results$r_value == final_results$r_mod_value) & (final_results$detected_value >= 0.5)
  final_results$tn = (final_results$r_value == final_results$r_mod_value) & (final_results$detected_value < 0.5)
  final_results$fn = (final_results$r_value != final_results$r_mod_value) & (final_results$detected_value < 0.5)
  
  #drop known deletions
  final_results <- subset(final_results, (final_results$chr != '10') & (final_results$chr != '13') & (final_results$chr != '14') & (final_results$chr != '19'))
  
  # Calculate precision and sensitivity
  TP = sum(final_results$tp)
  FP = sum(final_results$fp)
  FN = sum(final_results$fn)
  TN = sum(final_results$tn)
  precision <- TP / (TP + FP)
  sensitivity <- TP / (TP + FN)
  
  #print((precision))
  #print((sensitivity))
  return(c(sensitivity = sensitivity, precision = precision, tp=TP, fp=FP, fn=FN, tn=TN))
}

#
process_standard <- function(hb, r, r_simulated) {
  hb$setAlleleMats(r.init=r_simulated, n.sc.init=cov.sc, het.deviance.threshold=0.1, n.cores=1)
  hb$setGeneFactors(txdb) ## map SNPs to genes
  hb$calcAlleleCnvBoundariesSimple(init=TRUE, verbose=FALSE) ## HMM
  hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
  ## look at final results
  results <- hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
  performance <- process_results(r, r_simulated, results)
  return(performance)
}

process_PCA <- function(hb, r, r_simulated) {
  hb$setAlleleMats(r.init=r_simulated, n.sc.init=cov.sc, het.deviance.threshold=0.1, n.cores=1)
  hb$setGeneFactors(txdb) ## map SNPs to genes
  hb$calcAlleleCnvBoundariesPCA(init=TRUE, verbose=FALSE) ## HMM
  hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
  ## look at final results
  results <- hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
  performance <- process_results(r, r_simulated, results)
  return(performance)
}

process_DBSCAN <- function(hb, r, r_simulated) {
  hb$setAlleleMats(r.init=r_simulated, n.sc.init=cov.sc, het.deviance.threshold=0.1, n.cores=1)
  hb$setGeneFactors(txdb) ## map SNPs to genes
  hb$calcAlleleCnvBoundariesDBSCAN(init=TRUE, verbose=FALSE) ## HMM
  hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
  ## look at final results
  results <- hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
  performance <- process_results(r, r_simulated, results)
  return(performance)
}

process_old <- function(hb, r, r_simulated) {
  hb$setAlleleMats(r.init=r_simulated, n.sc.init=cov.sc, het.deviance.threshold=0.1, n.cores=1)
  hb$setGeneFactors(txdb) ## map SNPs to genes
  hb$calcAlleleCnvBoundariesOld(init=TRUE, verbose=FALSE) ## HMM
  hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
  ## look at final results
  results <- hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
  performance <- process_results(r, r_simulated, results)
  return(performance)
}

###
#Simulating deletes for cluster techniques
###

# List of known deletion sizes -> 40 translates to ~10Mb
del_sizes <- c(40, 80, 120, 160, 200, 240)

# Initialize a list to store the results
results <- list()

for (del_size in del_sizes) {
  cat("Processing deletion size:", del_size, "\n")
  
  # Initialize an empty data frame to store the results for this deletion size
  del_results <- data.frame(del_size = integer(), model = character(), sensitivity = numeric(), precision = numeric(),
                            tp = integer(), fp = integer(), fn = integer(), tn = integer(), stringsAsFactors = FALSE)
  
  
  hb <- new('HoneyBADGER', name='MGH31')
  data(r) ## alternate allele
  data(cov.sc) ## total coverage
  library(TxDb.Hsapiens.UCSC.hg19.knownGene) ## in order to map SNPs to genes
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  r_simulated <- simulate_deletions(r, del_size, 0.7)

  print("Running standard")
  result_standard <- tryCatch(process_standard(hb, r, r_simulated), error = function(e) {
    cat("Error in", "Standard", "\n")
    return(c(sensitivity = 0, precision = 0, tp=0, fp=0, fn=0, tn=0))
  })
  # Add the model results to the data frame
  del_results <- rbind(del_results, c(del_size, "Standard", result_standard))
  write.csv(del_results, file = paste0("tmp_sim_data_cluster/results_del_size_", del_size, "_v1.csv"), row.names = FALSE)
  
  print("Running PCA")
  result_PCA <- tryCatch(process_PCA(hb, r, r_simulated), error = function(e) {
    cat("Error in", "PCA", "\n")
    return(c(sensitivity = 0, precision = 0, tp=0, fp=0, fn=0, tn=0))
  })
  # Add the model results to the data frame
  del_results <- rbind(del_results, c(del_size, "PCA", result_PCA))
  write.csv(del_results, file = paste0("tmp_sim_data_cluster/results_del_size_", del_size, "_v2.csv"), row.names = FALSE)
  
  print("Running DBSCAN")
  result_DBSCAN <- tryCatch(process_DBSCAN(hb, r, r_simulated), error = function(e) {
    cat("Error in", "DBSCAN", "\n")
    return(c(sensitivity = 0, precision = 0, tp=0, fp=0, fn=0, tn=0))
  })
  # Add the model results to the data frame
  del_results <- rbind(del_results, c(del_size, "DBSCAN", result_DBSCAN))
  
  # Add the deletion size results to the list
  results[[as.character(del_size)]] <- del_results
  
  print("Finishing up deletion size - the results are below")
  print(results)
  
  # Save the deletion size results to a CSV file
  write.csv(del_results, file = paste0("tmp_sim_data_cluster/results_del_size_", del_size, "_v3.csv"), row.names = FALSE)
  
  print(result_standard)
  print(result_PCA)
  print(result_DBSCAN)
  
}




###
#Simulating deletes for allele normalization comparison
###

# List of known deletion sizes
del_sizes <- c(40, 80, 120, 160, 200, 240)

# Initialize a list to store the results
results <- list()

for (del_size in del_sizes) {
  cat("Processing deletion size:", del_size, "\n")
  
  # Initialize an empty data frame to store the results for this deletion size
  del_results <- data.frame(del_size = integer(), model = character(), sensitivity = numeric(), precision = numeric(),
                            tp = integer(), fp = integer(), fn = integer(), tn = integer(), stringsAsFactors = FALSE)
  
  
  hb <- new('HoneyBADGER', name='MGH31')
  data(r) ## alternate allele
  data(cov.sc) ## total coverage
  library(TxDb.Hsapiens.UCSC.hg19.knownGene) ## in order to map SNPs to genes
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  r_simulated <- simulate_deletions(r, del_size, 0.7)
  
  print("Running standard")
  result_standard <- tryCatch(process_standard(hb, r, r_simulated), error = function(e) {
    cat("Error in", "Standard", "\n")
    return(c(sensitivity = 0, precision = 0, tp=0, fp=0, fn=0, tn=0))
  })
  # Add the model results to the data frame
  del_results <- rbind(del_results, c(del_size, "Standard", result_standard))
  write.csv(del_results, file = paste0("tmp_sim_data_allele_norm/results_del_size_", del_size, "_v1.csv"), row.names = FALSE)
  
  print("Running old")
  result_old <- tryCatch(process_old(hb, r, r_simulated), error = function(e) {
    cat("Error in", "old", "\n")
    return(c(sensitivity = 0, precision = 0, tp=0, fp=0, fn=0, tn=0))
  })
  # Add the model results to the data frame
  del_results <- rbind(del_results, c(del_size, "Old", result_old))
  write.csv(del_results, file = paste0("tmp_sim_data_allele_norm/results_del_size_", del_size, "_v2.csv"), row.names = FALSE)
  
  # Add the deletion size results to the list
  results[[as.character(del_size)]] <- del_results
  
  print("Finishing up deletion size - the results are below")
  print(results)
  

}




###
#Helper code for investigating simulated matrix 
###

#diff_matrix <- r != r_simulated

# Count the number of differences for each cell
#diff_count <- colSums(diff_matrix)

# Calculate the percentage of differences for each cell
#total_elements <- nrow(r) * ncol(r)
#diff_percentage <- (diff_count / nrow(r)) * 100

# Create a summary data frame
#summary_df <- data.frame(Cell = colnames(r),
#                         Differences = diff_count,
#                         Percentage = diff_percentage)

# Print the summary data frame
#print(summary_df)



