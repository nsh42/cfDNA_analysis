# Define directory path containing read length output files created using BAM_file_fragment_sizes.R
directory <- "~/rds/hpc-work/Read_length_outputfiles"

file_list <- list.files(directory, pattern = ".txt",
                          full.names = TRUE)
  
samples <- character(length = length(file_list))
mean_frag_size <- vector(length = length(file_list))

counter = 1
 
# Loop through files and record sample name and mean fragment size
for (file in file_list) {
    
  print(paste("Scanning file: ", file, sep = ""))
    
  file_name <- basename(file)
  sample <- gsub(".txt", "", file_name)
  sample <- paste("P", sample, sep = "")
  samples[counter] <- sample
    
  fragment_lengths <- scan(file)
  fragment_lengths <- fragment_lengths[which(fragment_lengths <= 500)]
    
  sample_mean_frag_size <- sum(fragment_lengths)/length(fragment_lengths)
  mean_frag_size[counter] <- sample_mean_frag_size

  counter = counter + 1
    
}

mean_frag_sizes_df <- data.frame("sample" = samples, "mean_frag_size" = mean_frag_size)

write.csv(mean_frag_sizes_df, file = "~/rds/hpc-work/output_Rplots/mean_frag_size_table.csv", row.names = F)





