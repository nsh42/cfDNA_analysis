library(mclust)

# Define directory path containing read length output files created using BAM_file_fragment_sizes.R
directory <- "~/rds/hpc-work/Read_length_outputfiles"

file_list <- list.files(directory, pattern = ".txt",
                        full.names = TRUE)

for (file in file_list) {
  
  print(paste("Scanning file: ", file, sep = ""))
  
  file_name <- basename(file)
  sample <- gsub(".txt", "", file_name)
  sample <- paste("P", sample, sep = "")
  patient <- gsub("V.*", "", sample)
  
  fragment_lengths <- scan(file)
  fragment_lengths <- fragment_lengths[which(fragment_lengths <= 500)]
  
  # Unequal variance, 4 components
  model <- Mclust(fragment_lengths, G = 4, modelNames = "V")
  
  png(paste("~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/", patient, "/classification_plot_final_", sample,".png", sep = ""), width = 750, height = 500)
  plot(model, what = "classification", main = "Gaussian components in cfDNA modelled by Mclust")
  dev.off()
  
  density_G4 <- densityMclust(fragment_lengths, G = 4)
  
  png(paste("~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/", patient, "/density_plot_final_", sample,".png", sep = ""), width = 750, height = 500)
  plot(density_G4, data = fragment_lengths, what = "density", col = "blue", hist.col = "lightgrey", xlab = "Fragment length (bp)", 
  main = "Fit of density estimated by Mclust to histogram of observed data")
  
  dev.off()
  
}


