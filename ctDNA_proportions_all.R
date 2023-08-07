library(mclust)

# Define directory path containing read length output files created using BAM_file_fragment_sizes.R
directory <- "~/rds/hpc-work/Read_length_outputfiles"

file_list <- list.files(directory, pattern = ".txt",
                        full.names = TRUE)


ctDNA_proportions_df <- data.frame(samples = character(length(file_list)),
                                   mean1 = character(length(file_list)),
                                   mean2 = character(length(file_list)),
                                   mean3 = character(length(file_list)),
                                   mean4 = character(length(file_list)),
                                   prop1 = character(length(file_list)),
                                   prop2 = character(length(file_list)),
                                   prop3 = character(length(file_list)),
                                   prop4 = character(length(file_list)),
                                   var1 = character(length(file_list)),
                                   var2 = character(length(file_list)),
                                   var3 = character(length(file_list)),
                                   var4 = character(length(file_list)),
                                   ctDNA_proportions = numeric(length(file_list)))

counter = 1

for (file in file_list) {
  
  print(paste("Scanning file: ", file, sep = ""))
  
  file_name <- basename(file)
  sample <- gsub(".txt", "", file_name)
  sample <- paste("P", sample, sep = "")
  patient <- gsub("V.*", "", sample)
  
  fragment_lengths <- scan(file)
  fragment_lengths <- fragment_lengths[which(fragment_lengths <= 500)]
  
  # Unqual variance, 4 components
  model <- Mclust(fragment_lengths, G = 4, modelNames = "V")
  
  mean_values <- as.vector(model$parameters$mean)
  proportions <- as.vector(model$parameters$pro)
  variances <- as.vector(model$parameters$variance$sigmasq)
  
  if (mean_values[2] < 160){
    ctDNA_proportion <- proportions[2]
  } else if (mean_values[3] < 160){
    ctDNA_proportion <- proportions[3]
  } else {
    ctDNA_proportion <- 0
  }
  
  ctDNA_proportions_df[counter, "samples"] <- sample
  ctDNA_proportions_df[counter, "mean1"] <- mean_values[1]
  ctDNA_proportions_df[counter, "mean2"] <- mean_values[2]
  ctDNA_proportions_df[counter, "mean3"] <- mean_values[3]
  ctDNA_proportions_df[counter, "mean4"] <- mean_values[4]
  ctDNA_proportions_df[counter, "prop1"] <- proportions[1]
  ctDNA_proportions_df[counter, "prop2"] <- proportions[2]
  ctDNA_proportions_df[counter, "prop3"] <- proportions[3]
  ctDNA_proportions_df[counter, "prop4"] <- proportions[4]
  ctDNA_proportions_df[counter, "var1"] <- variances[1]
  ctDNA_proportions_df[counter, "var2"] <- variances[2]
  ctDNA_proportions_df[counter, "var3"] <- variances[3]
  ctDNA_proportions_df[counter, "var4"] <- variances[4]
  ctDNA_proportions_df[counter, "ctDNA_proportions"] <- ctDNA_proportion
  
  write.csv(ctDNA_proportions_df, file = "~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/ctDNA_proportions_all_data.csv", row.names = F)
  
  counter = counter + 1
}





