# ******** Plot fragment frequency for particular patients ************ #
directory <- "~/rds/hpc-work/Read_length_outputfiles"

library(tidyr)
library(ggplot2)

# Create function to plot fragment size distribution per patient
plot_patient_frag_sizes <- function(patient_number){
  
  ichorCNA_data <- read.csv("~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/ichorCNAOldFiles/TableS3.txt", sep = "\t")
  
  file_list <- list.files(directory, 
                          pattern = paste("^", patient_number, "V.*\\.txt$", sep = ""), 
                          full.names = TRUE)
  
  frag_data <- data.frame("Fragment_size" = 1:500)
  
  ichorCNA_values <- vector(length = length(file_list))
  Samples <- character(length = length(file_list))
  counter = 1
  
  for (file in file_list) {
    
    print(paste("Scanning file: ", file, sep = ""))
    
    file_name <- basename(file)
    sample <- gsub(".txt", "", file_name)
    sample <- paste("P", sample, sep = "")
    Samples[counter] <- sample
    
    fragment_lengths <- scan(file)
    fragment_lengths <- fragment_lengths[which(fragment_lengths <= 500)]
    
    fragment_freq <- table(fragment_lengths)
    
    fragment_freq_df <- data.frame("Fragment_size" = as.numeric(row.names(fragment_freq)))
    column_name <- sample
    fragment_freq_df[column_name] <- as.vector(fragment_freq)
    
    frag_data <- merge(frag_data, fragment_freq_df, by = "Fragment_size", all.x = TRUE)
    
    sample_ichorCNA <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "ichorCNA"]
    ichorCNA_values[counter] <- sample_ichorCNA
    counter = counter + 1
    
  }
  
  frag_data[is.na(frag_data)] <- 0
  
  for (i in 1:length(file_list)){
    column_sum <- sum(frag_data[,(i+1)])
    frag_data[,(i+1)]<-(frag_data[,(i+1)])/column_sum
  }
  
  frag_data_long <- tidyr::gather(frag_data, key = "Visit", value = "Frequency", -Fragment_size)
  frag_data_long$ichorCNA <- NA
  for (row in 1:nrow(frag_data_long)){
    Visit <- frag_data_long$Visit[row]
    frag_data_long$ichorCNA[row] <- ichorCNA_data[which(ichorCNA_data$SampleID == Visit), "ichorCNA"]
  }
  
  plot1 <- ggplot(frag_data_long, aes(x = Fragment_size, y = Frequency, color = ichorCNA, group = Visit)) +
    geom_line() +
    ggtitle(paste("Fragment size distributions for patient P", patient_number, sep = "")) +
    scale_color_gradient(low = "blue", high = "red", na.value = NA, guide = "colorbar") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 3, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title = element_text(size = 12),
          panel.grid = element_blank()) +
    xlab("Fragment Size") +
    ylab("Proportion of fragments")
  
  ggsave(filename = paste("~/rds/hpc-work/output_Rplots/Frag_size_dist_per_patient/P", patient_number, ".png", sep = ""), plot = plot1, height = 5, width = 9, device = "png")
  
  if (all(ichorCNA_values == 0) == FALSE){
    
    plot2 <- ggplot(frag_data_long, aes(x = Fragment_size, y = Frequency, color = ichorCNA, group = Visit)) +
      geom_line() +
      ggtitle(paste("Fragment size distributions for patient P", patient_number, sep = "")) +
      scale_color_gradient(low = "blue", high = "red", na.value = NA, guide = "legend", breaks = ichorCNA_values, labels = Samples) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 3, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.title = element_text(size = 12),
            panel.grid = element_blank()) +
      xlab("Fragment Size") +
      ylab("Proportion of fragments") +
      labs(color = "Patient visit")
    
    ggsave(filename = paste("~/rds/hpc-work/output_Rplots/Frag_size_dist_per_patient/Pervisit_P", patient_number, ".png", sep = ""), plot = plot2, height = 5, width = 9, device = "png")
    
  } else {
    print(paste("All ichorCNA values for patient P", patient_number, " are 0", sep = ""))
  }
  
}

all_files <- list.files(directory, 
                        pattern = "*\\.txt$", 
                        full.names = TRUE)

# Create vector containing all patient numbers
patient_numbers <- vector(length = length(all_files))
for (file in 1:length(all_files)){
  file_path = all_files[file]
  patient_number <- gsub("V.*\\.txt", "", basename(file_path))
  patient_number <- as.numeric(patient_number)
  patient_numbers[file] <- patient_number
}

patient_numbers <- unique(patient_numbers)

# Loop through patients and apply function to plot their fragment size distribution
for (i in 1:length(patient_numbers)){
  patient_number <- patient_numbers[i]
  plot_patient_frag_sizes(patient_number = patient_number)
}




