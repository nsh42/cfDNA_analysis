# ******** Plot binned fragment frequency histogram for all samples together ************ #

# Read in csv file containing frequency of fragments in each size bin, created using Find_binned_fragment_size_frequencies.R
fragment_frequencies <- as.data.frame(read.csv("~/rds/hpc-work/Read_length_outputfiles/Fragment_bins_10_frequency_table.csv", row.names = 1))
row.names(fragment_frequencies) <- paste("P", row.names(fragment_frequencies), sep = "")

bins <- seq(0, 500, 10)
colnames(fragment_frequencies) <- paste(bins[-length(bins)], bins[-1], sep = "-")

bin_sums <- colSums(fragment_frequencies)
bin_sum_df<-data.frame("Fragment.size.bin" = colnames(fragment_frequencies), 
                       "Frequency" = bin_sums)

bin_sum_df$Fragment.size.bin <- factor(bin_sum_df$Fragment.size.bin, levels = colnames(fragment_frequencies))

library("ggplot2")

total_frag_size_dist_plot <- ggplot(bin_sum_df) +
  geom_bar(aes(x = Fragment.size.bin, y = Frequency), stat = "identity", fill = "deepskyblue3") +
  ggtitle("Fragment size distribution across all samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid = element_blank()) +
  xlab("Fragment size (bp)") +
  ylab("Frequency")

ggsave(filename = "~/rds/hpc-work/output_Rplots/total_frag_size_dist10_plot.png", plot = total_frag_size_dist_plot, height = 5, width = 15, device = "png")


# ******** Plot binned fragment frequencies by ichorCNA range (binned by 10) ************ #

ichorCNA_data <- read.csv("~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/ichorCNAOldFiles/TableS3.txt", sep = "\t")

fragment_frequencies <- as.data.frame(read.csv("~/rds/hpc-work/Read_length_outputfiles/Fragment_bins_10_frequency_table.csv", row.names = 1))
row.names(fragment_frequencies) <- paste("P", row.names(fragment_frequencies), sep = "")

bins <- seq(0, 500, 10)
colnames(fragment_frequencies) <- paste(bins[-length(bins)], bins[-1], sep = "-")

bin_sums <- colSums(fragment_frequencies)
bin_sum_df<-data.frame("Fragment.size.bin" = colnames(fragment_frequencies), 
                       "Frequency" = bin_sums)

bin_sum_df$Fragment.size.bin <- factor(bin_sum_df$Fragment.size.bin, levels = colnames(fragment_frequencies))

fragment_frequencies$ichorCNA <- NA
for (row in 1:nrow(fragment_frequencies)){
  sample <- row.names(fragment_frequencies)[row]
  ichorCNA <- ichorCNA_data[which(ichorCNA_data$SampleID == sample), "ichorCNA"]
  if (ichorCNA <= 20){
    ichorCNA_range <- "0-20"
  } else if (ichorCNA > 20 && ichorCNA <= 40){
    ichorCNA_range <- "20-40"
  } else if (ichorCNA > 40 && ichorCNA <= 60){
    ichorCNA_range <- "40-60"
  } else {
    ichorCNA_range <- "60-80"
  }
  fragment_frequencies$ichorCNA[row] <- ichorCNA_range
}

fragment_frequencies$ichorCNA <- as.factor(fragment_frequencies$ichorCNA)

library(dplyr)
frequencies_by_ichorCNA <- fragment_frequencies %>%
  group_by(ichorCNA) %>%
  summarize(size_0_to_20 = sum(`0-10`),
            size_10_to_20 = sum(`10-20`),
            size_20_to_30 = sum(`20-30`),
            size_30_to_40 = sum(`30-40`),
            size_40_to_50 = sum(`40-50`),
            size_50_to_60 = sum(`50-60`),
            size_60_to_70 = sum(`60-70`),
            size_70_to_80 = sum(`70-80`),
            size_80_to_90 = sum(`80-90`),
            size_90_to_100 = sum(`90-100`),
            size_100_to_110 = sum(`100-110`),
            size_110_to_120 = sum(`110-120`),
            size_120_to_130 = sum(`120-130`),
            size_130_to_140 = sum(`130-140`),
            size_140_to_150 = sum(`140-150`),
            size_150_to_160 = sum(`150-160`),
            size_160_to_170 = sum(`160-170`),
            size_170_to_180 = sum(`170-180`),
            size_180_to_190 = sum(`180-190`),
            size_190_to_200 = sum(`190-200`),
            size_200_to_210 = sum(`200-210`),
            size_210_to_220 = sum(`210-220`),
            size_220_to_230 = sum(`220-230`),
            size_230_to_240 = sum(`230-240`),
            size_240_to_250 = sum(`240-250`),
            size_250_to_260 = sum(`250-260`),
            size_260_to_270 = sum(`260-270`),
            size_270_to_280 = sum(`270-280`),
            size_280_to_290 = sum(`280-290`),
            size_290_to_300 = sum(`290-300`),
            size_300_to_310 = sum(`300-310`),
            size_310_to_320 = sum(`310-320`),
            size_320_to_330 = sum(`320-330`),
            size_330_to_340 = sum(`330-340`),
            size_340_to_350 = sum(`340-350`),
            size_350_to_360 = sum(`350-360`),
            size_360_to_370 = sum(`360-370`),
            size_370_to_380 = sum(`370-380`),
            size_380_to_390 = sum(`380-390`),
            size_390_to_400 = sum(`390-400`),
            size_400_to_410 = sum(`400-410`),
            size_410_to_420 = sum(`410-420`),
            size_420_to_430 = sum(`420-430`),
            size_430_to_440 = sum(`430-440`),
            size_440_to_450 = sum(`440-450`),
            size_450_to_460 = sum(`450-460`),
            size_460_to_470 = sum(`460-470`),
            size_470_to_480 = sum(`470-480`),
            size_480_to_490 = sum(`480-490`),
            size_490_to_500 = sum(`490-500`))

row_sums <- rowSums(frequencies_by_ichorCNA[-1])
frequencies_by_ichorCNA<- as.data.frame(frequencies_by_ichorCNA)
frequencies_by_ichorCNA[1,-1] <- frequencies_by_ichorCNA[1,-1]/row_sums[1]
frequencies_by_ichorCNA[2,-1] <- frequencies_by_ichorCNA[2,-1]/row_sums[2]
frequencies_by_ichorCNA[3,-1] <- frequencies_by_ichorCNA[3,-1]/row_sums[3]
frequencies_by_ichorCNA[4,-1] <- frequencies_by_ichorCNA[4,-1]/row_sums[4]

library(tidyr)
melted_df <- tidyr::pivot_longer(frequencies_by_ichorCNA, cols = -ichorCNA, names_to = "Size_Bin", values_to = "Frequency")

melted_df$Size_Bin <- factor(melted_df$Size_Bin, levels = melted_df$Size_Bin[1:50])

frag_size_by_ichorCNA_plot <- ggplot(melted_df, aes(x = Size_Bin, y = Frequency, color = ichorCNA, group = ichorCNA), ) +
  geom_line() +
  ggtitle("Fragment size distribution by ichorCNA range") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"),
        legend.position = c(0.95,0.8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid = element_blank()) +
  xlab("Fragment size (bp)") +
  ylab("Proportion of fragments") +
  scale_x_discrete(labels = seq(0, 500, 10)) +
  scale_color_manual(values = c("blue", "purple", "mediumvioletred", "red")) 

ggsave(filename = "~/rds/hpc-work/output_Rplots/Frag_size_dist_by_ichorCNA.png", plot = frag_size_by_ichorCNA_plot, height = 5, width = 15, device = "png")


