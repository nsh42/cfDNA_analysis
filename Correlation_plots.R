# *********** SLFR vs ichorCNA correlation ************ #

ichorCNA_data <- read.csv("~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/ichorCNAOldFiles/TableS3.txt", sep = "\t")
fragment_frequencies <- as.data.frame(read.csv("~/rds/hpc-work/Read_length_outputfiles/Fragment_bins_10_frequency_table.csv", row.names = 1))

small_fragments <- rowSums(fragment_frequencies[,c(9:16)]) + fragment_frequencies[,8]*7/10 + fragment_frequencies[,17]*5/10  # define small fragments (73-165bp)
large_fragments <- rowSums(fragment_frequencies[,c(18:25)]) + fragment_frequencies[,17]*5/10 + fragment_frequencies[,26]*3/10 # define long fragments (166-253bp)

small_over_large <- small_fragments/large_fragments

fragsize_ichorCNA_df<-data.frame("small_over_large" = small_over_large)
row.names(fragsize_ichorCNA_df) <- paste("P", row.names(fragment_frequencies), sep = "")

fragsize_ichorCNA_df$ichorCNA <- NA
for (row in 1:nrow(fragsize_ichorCNA_df)){
  sample <- row.names(fragsize_ichorCNA_df)[row]
  ichorCNA <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "ichorCNA"]
  fragsize_ichorCNA_df$ichorCNA[row] <- ichorCNA
}

fit <- lm(fragsize_ichorCNA_df$small_over_large ~ fragsize_ichorCNA_df$ichorCNA)

fragsize_ichorCNA_corr_plot <- ggplot(fragsize_ichorCNA_df, aes(x = ichorCNA, y = small_over_large)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between SLFR and ichorCNA",
    x = "ichorCNA",
    y = "Ratio of small to large fragments"
  )

ggsave(filename = "~/rds/rds/hpc-work/output_Rplots/fragsize_ichorCNA_corr_plot.png",
       fragsize_ichorCNA_corr_plot, width = 9, height = 5, device = png)

cor.test(fragsize_ichorCNA_df$small_over_large, fragsize_ichorCNA_df$ichorCNA, method = "pearson") 
cor.test(fragsize_ichorCNA_df$small_over_large, fragsize_ichorCNA_df$ichorCNA, method = "spearman") 


# *********** SLFR vs entropy correlation ************ #
fragment_frequencies <- as.data.frame(read.csv("~/rds/hpc-work/Read_length_outputfiles/Fragment_bins_10_frequency_table.csv", row.names = 1))
entropy_df <- read.csv("~/rds/hpc-work/output_Rplots/entropy_per_sample.csv", sep = ",")

small_fragments <- rowSums(fragment_frequencies[,c(9:16)]) + fragment_frequencies[,8]*7/10 + fragment_frequencies[,17]*5/10 # define small fragments (73-165bp)
large_fragments <- rowSums(fragment_frequencies[,c(18:25)]) + fragment_frequencies[,17]*5/10 + fragment_frequencies[,26]*3/10 # define long fragments (166-253bp)

small_over_large <- small_fragments/large_fragments

sizeratio_entropy_df<-data.frame("small_over_large" = small_over_large)
row.names(sizeratio_entropy_df) <- paste("P", row.names(fragment_frequencies), sep = "")

sizeratio_entropy_df$entropy <- NA
for (row in 1:nrow(sizeratio_entropy_df)){
  sample <- row.names(sizeratio_entropy_df)[row]
  entropy <- entropy_df[which(entropy_df$sample == sample), "entropy"]
  sizeratio_entropy_df$entropy[row] <- entropy
}

fit <- lm(sizeratio_entropy_df$small_over_large ~ sizeratio_entropy_df$entropy)

fragsize_entropy_corr_plot <- ggplot(sizeratio_entropy_df, aes(x = entropy, y = small_over_large)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between SLFR and entropy",
    x = "Entropy",
    y = "Ratio of small to large fragments"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/fragsize_entropy_corr_plot.png",
       fragsize_entropy_corr_plot, width = 9, height = 5, device = png)

cor.test(sizeratio_entropy_df$small_over_large, sizeratio_entropy_df$entropy, method = "pearson") 
cor.test(sizeratio_entropy_df$small_over_large, sizeratio_entropy_df$entropy, method = "spearman") 


# *********** SLFR vs MFS correlation ************ #
fragment_frequencies <- as.data.frame(read.csv("~/rds/hpc-work/Read_length_outputfiles/Fragment_bins_10_frequency_table.csv", row.names = 1))
mean_fragment_sizes <- read.csv("~/rds/hpc-work/output_Rplots/mean_frag_size_table.csv", sep = ",")

small_fragments <- rowSums(fragment_frequencies[,c(9:16)]) + fragment_frequencies[,8]*7/10 + fragment_frequencies[,17]*5/10
large_fragments <- rowSums(fragment_frequencies[,c(18:25)]) + fragment_frequencies[,17]*5/10 + fragment_frequencies[,26]*3/10

small_over_large <- small_fragments/large_fragments

sizeratio_meansize_df<-data.frame("small_over_large" = small_over_large)
row.names(sizeratio_meansize_df) <- paste("P", row.names(fragment_frequencies), sep = "")

sizeratio_meansize_df$mean_size <- NA
for (row in 1:nrow(sizeratio_meansize_df)){
  sample <- row.names(sizeratio_meansize_df)[row]
  mean_size <- mean_fragment_sizes[which(mean_fragment_sizes$sample == sample), "mean_frag_size"]
  sizeratio_meansize_df$mean_size[row] <- mean_size
}

fit <- lm(sizeratio_meansize_df$small_over_large ~ sizeratio_meansize_df$mean_size)

fragsize_meansize_corr_plot <- ggplot(sizeratio_meansize_df, aes(x = mean_size, y = small_over_large)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between SLFR and MFS",
    x = "mean fragment size (bp)",
    y = "Ratio of small to large fragments"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/fragsize_meansize_corr_plot.png",
       fragsize_meansize_corr_plot, width = 9, height = 5, device = png)

cor.test(sizeratio_meansize_df$small_over_large, sizeratio_meansize_df$mean_size, method = "pearson") 
cor.test(sizeratio_meansize_df$small_over_large, sizeratio_meansize_df$mean_size, method = "spearman") 


# ********** MFS vs ichorCNA correlation ************** #
mean_fragment_sizes <- read.csv("~/rds/hpc-work/output_Rplots/mean_frag_size_table.csv", sep = ",")
ichorCNA_data <- read.csv("~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/ichorCNAOldFiles/TableS3.txt", sep = "\t")

mean_fragment_sizes$ichorCNA <- NA
for (i in 1:nrow(mean_fragment_sizes)){
  sample <- mean_fragment_sizes$sample[i]
  ichorCNA <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "ichorCNA"]
  mean_fragment_sizes$ichorCNA[i] <- ichorCNA
}

fit <- lm(mean_fragment_sizes$mean_frag_size ~ mean_fragment_sizes$ichorCNA)

mean_fragsize_ichorCNA_corr_plot <- ggplot(mean_fragment_sizes, aes(x = ichorCNA, y = mean_frag_size)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between MFS and ichorCNA",
    x = "ichorCNA",
    y = "Mean fragment size (bp)"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/mean_fragsize_ichorCNA_corr_plot.png",
       mean_fragsize_ichorCNA_corr_plot, width = 9, height = 5, device = png)

cor.test(mean_fragment_sizes$mean_frag_size, mean_fragment_sizes$ichorCNA, method = "pearson") 
cor.test(mean_fragment_sizes$mean_frag_size, mean_fragment_sizes$ichorCNA, method = "spearman") 


# ********** entropy vs ichorCNA correlation ************** #
entropy_df <- read.csv("~/rds/hpc-work/output_Rplots/entropy_per_sample.csv", sep = ",")
ichorCNA_data <- read.csv("~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/ichorCNAOldFiles/TableS3.txt", sep = "\t")

entropy_df$ichorCNA <- NA
for (i in 1:nrow(entropy_df)){
  sample <- entropy_df$sample[i]
  ichorCNA <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "ichorCNA"]
  entropy_df$ichorCNA[i] <- ichorCNA
}

fit <- lm(entropy_df$entropy ~ entropy_df$ichorCNA)

entropy_ichorCNA_corr_plot <- ggplot(entropy_df, aes(x = ichorCNA, y = entropy)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between entropy and ichorCNA",
    x = "ichorCNA",
    y = "Entropy"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/entropy_ichorCNA_corr_plot_notbinned.png",
       mean_fragsize_ichorCNA_corr_plot, width = 9, height = 5, device = png)

cor.test(entropy_df$entropy, entropy_df$ichorCNA, method = "pearson")
cor.test(entropy_df$entropy, entropy_df$ichorCNA, method = "spearman") 



# ********** MFS vs entropy correlation ************** #
entropy_df <- read.csv("~/rds/hpc-work/output_Rplots/entropy_per_sample.csv", sep = ",")
mean_fragment_sizes <- read.csv("~/rds/hpc-work/output_Rplots/mean_frag_size_table.csv", sep = ",")

entropy_df$mean_size <- NA
for (i in 1:nrow(entropy_df)){
  sample <- entropy_df$sample[i]
  mean_size <- mean_fragment_sizes[which(mean_fragment_sizes$sample == sample), "mean_frag_size"]
  entropy_df$mean_size[i] <- mean_size
}

fit <- lm(entropy_df$entropy ~ entropy_df$mean_size)

entropy_meansize_corr_plot <- ggplot(entropy_df, aes(x = mean_size, y = entropy)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between entropy and mean fragment size",
    x = "Mean fragment size (bp)",
    y = "Entropy"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/entropy_meansize_corr_plot.png",
       entropy_meansize_corr_plot, width = 9, height = 5, device = png)

cor.test(entropy_df$entropy, entropy_df$mean_size, method = "pearson") 
cor.test(entropy_df$entropy, entropy_df$mean_size, method = "spearman") 


# *********** mclust ctDNA proportion vs ichorCNA correlation ************ #

ctDNA_proportion_df <- as.data.frame(read.csv("~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/ctDNA_proportions_all_data.csv"))

ctDNA_proportion_df$ichorCNA <- NA
for (row in 1:nrow(ctDNA_proportion_df)){
  sample <- ctDNA_proportion_df$samples[row]
  ichorCNA <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "ichorCNA"]
  ctDNA_proportion_df$ichorCNA[row] <- ichorCNA
  if(ctDNA_proportion_df$mean2[row] > 160 & ctDNA_proportion_df$mean2[row] < 165){
    ctDNA_proportion_df$ctDNA_proportions[row] <- ctDNA_proportion_df$prop2[row]
  }
}

fit <- lm(ctDNA_proportion_df$ctDNA_proportions ~ ctDNA_proportion_df$ichorCNA)

ctDNAprop_ichorCNA_corr_plot <- ggplot(ctDNA_proportion_df, aes(x = ichorCNA, y = ctDNA_proportions)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between ctDNA proportion and ichorCNA",
    x = "ichorCNA",
    y = "ctDNA proportion"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/ctDNAprop_ichorCNA_corr_plot_165.png",
       ctDNAprop_ichorCNA_corr_plot, width = 9, height = 5, device = png)

cor.test(ctDNA_proportion_df$ctDNA_proportions, ctDNA_proportion_df$ichorCNA, method = "pearson") 
cor.test(ctDNA_proportion_df$ctDNA_proportions, ctDNA_proportion_df$ichorCNA, method = "spearman") 


# *********** mclust ctDNA proportion vs MFS correlation ************ #

ctDNA_proportion_df <- as.data.frame(read.csv("~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/ctDNA_proportions_all_data.csv"))
mean_fragment_sizes <- read.csv("~/rds/hpc-work/output_Rplots/mean_frag_size_table.csv", sep = ",")

ctDNA_proportion_df$mean_frag_size <- NA
for (row in 1:nrow(ctDNA_proportion_df)){
  sample <- ctDNA_proportion_df$samples[row]
  mean_frag_size <- mean_fragment_sizes[which(mean_fragment_sizes$sample == sample), "mean_frag_size"]
  ctDNA_proportion_df$mean_frag_size[row] <- mean_frag_size
  if(ctDNA_proportion_df$mean2[row] > 160 & ctDNA_proportion_df$mean2[row] < 165){
    ctDNA_proportion_df$ctDNA_proportions[row] <- ctDNA_proportion_df$prop2[row]
  }
}

fit <- lm(ctDNA_proportion_df$ctDNA_proportions ~ ctDNA_proportion_df$mean_frag_size)

ctDNAprop_mean_frag_size_corr_plot <- ggplot(ctDNA_proportion_df, aes(x = mean_frag_size, y = ctDNA_proportions)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between ctDNA proportion and mean fragment size",
    x = "mean fragment size (bp)",
    y = "ctDNA proportion"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/ctDNAprop_meanfragsize_corr_plot_165.png",
       ctDNAprop_mean_frag_size_corr_plot, width = 9, height = 5, device = png)

cor.test(ctDNA_proportion_df$ctDNA_proportions, ctDNA_proportion_df$mean_frag_size, method = "pearson") 
cor.test(ctDNA_proportion_df$ctDNA_proportions, ctDNA_proportion_df$mean_frag_size, method = "spearman") 


# *********** mclust ctDNA proportion vs entropy correlation ************ #

ctDNA_proportion_df <- as.data.frame(read.csv("~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/ctDNA_proportions_all_data.csv"))
entropy_df <- read.csv("~/rds/hpc-work/output_Rplots/entropy_per_sample.csv", sep = ",")

ctDNA_proportion_df$entropy <- NA
for (row in 1:nrow(ctDNA_proportion_df)){
  sample <- ctDNA_proportion_df$samples[row]
  entropy <- entropy_df[which(entropy_df$sample == sample), "entropy"]
  ctDNA_proportion_df$entropy[row] <- entropy
  if(ctDNA_proportion_df$mean2[row] > 160 & ctDNA_proportion_df$mean2[row] < 165){
    ctDNA_proportion_df$ctDNA_proportions[row] <- ctDNA_proportion_df$prop2[row]
  }
}

fit <- lm(ctDNA_proportion_df$ctDNA_proportions ~ ctDNA_proportion_df$entropy)

ctDNAprop_entropy_corr_plot <- ggplot(ctDNA_proportion_df, aes(x = entropy, y = ctDNA_proportions)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between ctDNA proportion and entropy",
    x = "entropy",
    y = "ctDNA proportion"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/ctDNAprop_entropy_corr_plot_165.png",
       ctDNAprop_entropy_corr_plot, width = 9, height = 5, device = png)

cor.test(ctDNA_proportion_df$ctDNA_proportions, ctDNA_proportion_df$entropy, method = "pearson") 
cor.test(ctDNA_proportion_df$ctDNA_proportions, ctDNA_proportion_df$entropy, method = "spearman")



# *********** mclust ctDNA proportion vs small to large ratio correlation ************ #

ctDNA_proportion_df <- as.data.frame(read.csv("~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/ctDNA_proportions_all_data.csv"))
fragment_frequencies <- as.data.frame(read.csv("~/rds/hpc-work/Read_length_outputfiles/Fragment_bins_10_frequency_table.csv", row.names = 1))

small_fragments <- rowSums(fragment_frequencies[,c(9:16)]) + fragment_frequencies[,8]*7/10 + fragment_frequencies[,17]*5/10
large_fragments <- rowSums(fragment_frequencies[,c(18:25)]) + fragment_frequencies[,17]*5/10 + fragment_frequencies[,26]*3/10

small_over_large <- small_fragments/large_fragments
fragsize_ctDNA_prop_df<-data.frame("small_over_large" = small_over_large)
row.names(fragsize_ctDNA_prop_df) <- paste("P", row.names(fragment_frequencies), sep = "")

ctDNA_proportion_df$fragsize_ratio <- NA
for (row in 1:nrow(ctDNA_proportion_df)){
  sample <- ctDNA_proportion_df$samples[row]
  fragsize_ratio <- fragsize_ctDNA_prop_df[which(row.names(fragsize_ctDNA_prop_df) == sample), "small_over_large"]
  ctDNA_proportion_df$fragsize_ratio[row] <- fragsize_ratio
  if(ctDNA_proportion_df$mean2[row] > 160 & ctDNA_proportion_df$mean2[row] < 165){
    ctDNA_proportion_df$ctDNA_proportions[row] <- ctDNA_proportion_df$prop2[row]
  }
}

fit <- lm(ctDNA_proportion_df$ctDNA_proportions ~ ctDNA_proportion_df$fragsize_ratio)

ctDNAprop_fragsizeratio_corr_plot <- ggplot(ctDNA_proportion_df, aes(x = fragsize_ratio, y = ctDNA_proportions)) +
  geom_point(size = 0.5) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "red") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 12)) +
  labs(
    title = "Correlation between ctDNA proportion and ratio of small to large fragments",
    x = "Ratio of small to large fragments",
    y = "ctDNA proportion"
  )

ggsave(filename = "~/rds/hpc-work/output_Rplots/ctDNAprop_fragsizeratio_corr_plot_165.png",
       ctDNAprop_fragsizeratio_corr_plot, width = 9, height = 5, device = png)

cor.test(ctDNA_proportion_df$ctDNA_proportions, ctDNA_proportion_df$fragsize_ratio, method = "pearson") 
cor.test(ctDNA_proportion_df$ctDNA_proportions, ctDNA_proportion_df$fragsize_ratio, method = "spearman") 

