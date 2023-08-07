plot_gaussian_curve <- function(sample, patient, n, mean_values, variance_values, proportions, colors){
  
  # Create a function for the Gaussian probability density function
  gaussian <- function(x, mean, variance) {
    exp(-((x - mean)^2) / (2 * variance)) / sqrt(2 * pi * variance)
  }
  
  # Generate x-axis values for plotting
  x <- seq(0, 500, 10)
  
  data <- data.frame(x = rep(x, n),
                     y = unlist(lapply(1:n, function(i) {
                       proportions[i] * dnorm(x, mean = mean_values[i], sd = sqrt(variance_values[i]))
                     })),
                     curve = rep(1:n, each = length(x)))
  
  # Plot the Gaussian curves using ggplot
  plot <- ggplot(data, aes(x = x, y = y, color = factor(curve))) +
    geom_line(size = 1) +
    scale_color_manual(values = colors) +
    labs(x = "Fragment length (bp)", y = "Probability Density", title = "Mclust Gaussian Curves") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title = element_text(size = 12),
          legend.position = "none") 
  
  ggsave(filename = paste("~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/", patient, "/Recoloured_", sample, "_gaussian_curves.png", sep = ""),
         plot, width = 9, height = 5, device = png)
  
}


library(mclust)
library(ggplot2)

# Define directory path containing read length output files created using BAM_file_fragment_sizes.R
directory <- "~/rds/hpc-work/Read_length_outputfiles"

file_list <- list.files(directory, pattern = ".txt",
                        full.names = TRUE)

for (file in file_list[c(30:32, 63:66, 111:114, 131:134, 140:148, 157:160)]) {
  
  print(paste("Scanning file: ", file, sep = ""))
  
  file_name <- basename(file)
  sample <- gsub(".txt", "", file_name)
  sample <- paste("P", sample, sep = "")
  patient <- gsub("V.*", "", sample)
  
  fragment_lengths <- scan(file)
  fragment_lengths <- fragment_lengths[which(fragment_lengths <= 500)]
  
  # Unqual variance, 4 components
  model <- Mclust(fragment_lengths, G = 4, modelNames = "V")
  
  components <- unique(model$classification)
  n <- length(components)
  mean_values <- as.vector(model$parameters$mean)
  variance_values <- model$parameters$variance$sigmasq[components]
  proportions <- model$parameters$pro[components]
  
  colors <- character(length = 4)
  colors[1] <- "steelblue2"
  colors[4] <- "mediumpurple"
  if (mean_values[2] < 158){
      colors[2] <- "red"
  } else {
      colors[2] <- "plum1"
  }
  
  if (mean_values[3] < 158){
    colors[3] <- "red"
  } else {
    colors[3] <- "green"
  }
  
  mean_values <- mean_values[components]
  colors <- colors[components]
  
  plot_gaussian_curve(sample, patient, n, mean_values, variance_values, proportions, colors)
  
}

