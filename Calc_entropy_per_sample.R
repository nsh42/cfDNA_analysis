# ******** Calculate entropy per sample ************ #

# Define directory path containing read length output files created using BAM_file_fragment_sizes.R
directory <- "~/rds/hpc-work/Read_length_outputfiles"
  
file_list <- list.files(directory, 
                          pattern = ".txt", 
                          full.names = TRUE)
  
samples <- character(length = length(file_list))
entropy <- vector(length = length(file_list))
counter = 1

# Loop through files and create dataframe with sample name and entropy
for (file in file_list) {
    
    print(paste("Scanning file: ", file, sep = ""))
    
    file_name <- basename(file)
    sample <- gsub(".txt", "", file_name)
    sample <- paste("P", sample, sep = "")
    samples[counter] <- sample
    
    fragment_lengths <- scan(file)
    fragment_lengths <- fragment_lengths[which(fragment_lengths <= 500)]
    
    fragment_freq <- table(fragment_lengths)
    fragment_freq_df <- as.data.frame(fragment_freq)
    
    N <- length(fragment_lengths)
    p_is <- fragment_freq_df$Freq / N
    log_p_is <- log(p_is)
    p_is_log_p_is <- p_is * log_p_is
    H_B <- -1 * sum(p_is_log_p_is)
    
    entropy[counter] <- H_B
    counter = counter + 1
}

  
entropy_df <- data.frame("sample" = samples, "entropy" = entropy)

write.csv(entropy_df, file = "~/rds/hpc-work/output_Rplots/entropy_per_sample.csv", row.names = F)







