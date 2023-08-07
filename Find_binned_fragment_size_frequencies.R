# Define directory path containing fragment length txt files created using BAM_file_fragment_sizes.R
directory <- "~/rds/hpc-work/Read_length_outputfiles"

# Create list of all txt files in the directory
files <- list.files(directory, pattern = "*.txt", full.names = TRUE)

# Define fragment size bins
bins <- seq(0, 500, 10)

# Create an empty table to store the results
results <- matrix(0, nrow = length(files), ncol = length(bins) - 1)
colnames(results) <- paste(bins[-length(bins)], bins[-1], sep = "-")
rownames(results) <- tools::file_path_sans_ext(basename(files))

# Loop through each file
for (i in seq_along(files)) {
  
  print(paste("Scanning file: ", files[i], sep = ""))
  
  # Read the fragment sizes from the file
  fragment_sizes <- scan(files[i], quiet = TRUE)
  
  # Count the frequency of fragments in each size bin
  frequencies <- table(cut(fragment_sizes, breaks = bins, right = FALSE))
  
  # Store the frequencies in the results table
  results[i, ] <- frequencies
}

# Write the results to a csv file
write.csv(results, paste(directory, "/Fragment_bins_10_frequency_table.csv", sep = ""))
