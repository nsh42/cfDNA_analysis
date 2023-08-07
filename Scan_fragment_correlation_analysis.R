# Read in all data
Scan_prog_data <- read.csv("~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/ichorCNAOldFiles/TableS2.txt", sep = "\t")
mean_fragment_sizes <- read.csv("~/rds/hpc-work/output_Rplots/mean_frag_size_table.csv", sep = ",")
ichorCNA_data <- read.csv("~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/ichorCNAOldFiles/TableS3.txt", sep = "\t")
entropy_df <- read.csv("~/rds/hpc-work/output_Rplots/entropy_per_sample.csv", sep = ",")
ctDNA_proportion_df <- as.data.frame(read.csv("~/rds/hpc-work/output_Rplots/Mclust_plots/exploratory_analysis/ctDNA_proportions_all_data.csv"))

for (row in 1:nrow(ctDNA_proportion_df)){
  if(ctDNA_proportion_df$mean2[row] > 160 & ctDNA_proportion_df$mean2[row] < 165){
    ctDNA_proportion_df$ctDNA_proportions[row] <- ctDNA_proportion_df$prop2[row]
  }
}

fragment_frequencies <- as.data.frame(read.csv("~/rds/hpc-work/Read_length_outputfiles/Fragment_bins_10_frequency_table.csv", row.names = 1))
small_fragments <- rowSums(fragment_frequencies[,c(11:15)])
large_fragments <- rowSums(fragment_frequencies[,c(16:22)])
small_over_large <- small_fragments/large_fragments


Day <- vector(length = nrow(mean_fragment_sizes))
Patient <- vector(length = nrow(mean_fragment_sizes))
IchorCNA <- vector(length = nrow(mean_fragment_sizes))
for (row in 1:nrow(mean_fragment_sizes)){
  sample <- mean_fragment_sizes$sample[row]
  day <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "Day"]
  patient <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "Public.ID"]
  ichorCNA <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "ichorCNA"]
  Day[row] <- day
  Patient[row] <- patient
  IchorCNA[row] <- ichorCNA
}


Scan_dist_metrics_df <- data.frame("Patient" = Patient, 
                               "Sample" = mean_fragment_sizes$sample,
                               "Day" = Day,
                               "ichorCNA" = IchorCNA,
                               "Mean_frag_size" = mean_fragment_sizes$mean_frag_size,
                               "Small_to_large" = small_over_large,
                               "Entropy" = entropy_df$entropy,
                               "ctDNA_proportion" = ctDNA_proportion_df$ctDNA_proportions)

Patients <- unique(Scan_dist_metrics_df$Patient)

Scan_dist_metrics_df$ScanDay <- NA
Scan_dist_metrics_df$Progression <- NA
for (i in 1:length(Patients)){
  patient <- Patients[i]
  rows <- which(Scan_dist_metrics_df$Patient == patient)
  visit_days <- Scan_dist_metrics_df[which(Scan_dist_metrics_df$Patient == patient),"Day"]
  scan_days <- Scan_prog_data[which(Scan_prog_data$Public.ID == patient),"Day"]
  closest_scan <- sapply(visit_days, function(x) scan_days[which.min(abs(x - scan_days))])
  Scan_dist_metrics_df$ScanDay[rows] <- closest_scan
  Scan_prog <- character(length = length(closest_scan))
  for (i in 1:length(closest_scan)){
    scan <- closest_scan[i]
    scan_prog <- Scan_prog_data[which(Scan_prog_data$Day == scan & Scan_prog_data$Public.ID == patient),"Progression"]
    Scan_prog[i] <- scan_prog
  }
  Scan_dist_metrics_df$Progression[rows] <- Scan_prog
}



for (row in 1:nrow(Scan_dist_metrics_df)){
  visit_day <- Scan_dist_metrics_df$Day[row]
  scan_day <- Scan_dist_metrics_df$ScanDay[row]
  diff <- abs(visit_day - scan_day)
  if (diff > 90){
    Scan_dist_metrics_df$ScanDay[row] <- NA
    Scan_dist_metrics_df$Progression[row] <- NA
  }
}

Scan_dist_metrics_df <- Scan_dist_metrics_df[-c(which(is.na(Scan_dist_metrics_df$ScanDay))), ]

# IchorCNA #
plot1 <- ggplot(Scan_dist_metrics_df, aes(x = Progression, y = ichorCNA, fill = Progression)) +
  geom_boxplot() +
  labs(x = "Disease progression", y = "ichorCNA") +
  scale_fill_manual(values = c("seagreen3", "orangered")) + 
  theme_bw() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"), 
        axis.title = element_text(size = 12)) +
  ylim(c(0, 80)) +
  ggtitle("Correlation between disease progression and ichorCNA")

ggsave(filename = "~/rds/hpc-work/output_Rplots/Scan_correlation_plots/Scan_ichorCNA_boxplot.png",
       plot1, width = 9, height = 5, device = png)


# Mean fragment size #
plot2 <- ggplot(Scan_dist_metrics_df, aes(x = Progression, y = Mean_frag_size, fill = Progression)) +
  geom_boxplot() +
  labs(x = "Disease progression", y = "Mean fragment size (bp)") +
  scale_fill_manual(values = c("seagreen3", "orangered")) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"), 
        axis.title = element_text(size = 12)) +
  ylim(c(150, 200)) +
  ggtitle("Correlation between disease progression and mean fragment size")

ggsave(filename = "~/rds/hpc-work/output_Rplots/Scan_correlation_plots/Scan_meanfragsize_boxplot.png",
       plot2, width = 9, height = 5, device = png)


# Entropy #
plot3 <- ggplot(Scan_dist_metrics_df, aes(x = Progression, y = Entropy, fill = Progression)) +
  geom_boxplot() +
  labs(x = "Disease progression", y = "Entropy") +
  scale_fill_manual(values = c("seagreen3", "orangered")) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"), 
        axis.title = element_text(size = 12)) +
  ylim(c(4.4, 5.4)) +
  ggtitle("Correlation between disease progression and entropy")

ggsave(filename = "~/rds/hpc-work/output_Rplots/Scan_correlation_plots/Scan_entropy_boxplot.png",
       plot3, width = 9, height = 5, device = png)


# Small to large fragment ratio #
plot4 <- ggplot(Scan_dist_metrics_df, aes(x = Progression, y = Small_to_large, fill = Progression)) +
  geom_boxplot() +
  labs(x = "Disease progression", y = "small to large fragment ratio") +
  scale_fill_manual(values = c("seagreen3", "orangered")) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"), 
        axis.title = element_text(size = 12)) +
  ylim(c(0.1, 0.7)) +
  ggtitle("Correlation between disease progression and small to large fragment ratio")

ggsave(filename = "~/rds/hpc-work/output_Rplots/Scan_correlation_plots/Scan_smalltolarge_boxplot.png",
       plot4, width = 9, height = 5, device = png)


# ctDNA proportion #
plot5 <- ggplot(Scan_dist_metrics_df, aes(x = Progression, y = ctDNA_proportion, fill = Progression)) +
  geom_boxplot() +
  labs(x = "Disease progression", y = "Mclust ctDNA proportion") +
  scale_fill_manual(values = c("seagreen3", "orangered")) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"), 
        axis.title = element_text(size = 12)) +
  ylim(c(0, 0.35)) +
  ggtitle("Correlation between disease progression and GMM-predicted ctDNA proportion")

ggsave(filename = "~/rds/hpc-work/output_Rplots/Scan_correlation_plots/Scan_ctDNAprop_boxplot.png",
       plot5, width = 9, height = 5, device = png)


kruskal.test(Mean_frag_size ~ Progression, data = Scan_dist_metrics_df) 
kruskal.test(Small_to_large ~ Progression, data = Scan_dist_metrics_df) 
kruskal.test(Entropy ~ Progression, data = Scan_dist_metrics_df) 
kruskal.test(ichorCNA ~ Progression, data = Scan_dist_metrics_df) 
kruskal.test(ctDNA_proportion ~ Progression, data = Scan_dist_metrics_df)



# Logistic regression
install.packages("lmtest")
library("lmtest")

Scan_dist_metrics_df$Progression <- factor(Scan_dist_metrics_df$Progression, levels = c("YES", "NO"))

full_model <- glm(Progression ~ ichorCNA + Small_to_large + ctDNA_proportion + Mean_frag_size + Entropy, data = Scan_dist_metrics_df, family = binomial)

reduced_model <- glm(Progression ~ ichorCNA, data = Scan_dist_metrics_df, family = binomial)

lr_test <- lrtest(full_model, reduced_model)
print(lr_test)


