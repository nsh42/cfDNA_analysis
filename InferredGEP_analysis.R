library(ggplot2)


# ****** Investigating NDR coverage : merged files ****** #
NDR_cov <- read.csv("~/rds/hpc-work/EPICseq_all_merged_bams_output/coverage.ndr.matrix.merged.by.10.txt", sep = "\t")
NDR_cov <- NDR_cov[, grep("NDR_all", colnames(NDR_cov))]

NDR_cov_vector<-as.vector(unlist(NDR_cov))
median(NDR_cov_vector) 
range(NDR_cov_vector) 


# ****** Investigating 2k coverage : merged files ****** #

Cov_2k_gr0 <- read.csv("~/rds/hpc-work/EPICseq_all_merged_bams_output/coverage.2k.matrix.merged.by.0.txt", sep = "\t")
Cov_2k_gr0 <- Cov_2k_gr0[, grep("NDR_all", colnames(Cov_2k_gr0))]
Cov_2k_vector_gr0<-as.vector(unlist(Cov_2k_gr0))

Cov_2k_gr20 <- read.csv("~/rds/hpc-work/EPICseq_all_merged_bams_output/coverage.2k.matrix.merged.by.20.txt", sep = "\t")
Cov_2k_gr20 <- Cov_2k_gr20[, grep("NDR_all", colnames(Cov_2k_gr20))]
Cov_2k_vector_gr20<-as.vector(unlist(Cov_2k_gr20))

Cov_2k_gr40 <- read.csv("~/rds/hpc-work/EPICseq_all_merged_bams_output/coverage.2k.matrix.merged.by.40.txt", sep = "\t")
Cov_2k_gr40 <- Cov_2k_gr40[, grep("NDR_all", colnames(Cov_2k_gr40))]
Cov_2k_vector_gr40<-as.vector(unlist(Cov_2k_gr40))

Cov_2k_gr60 <- read.csv("~/rds/hpc-work/EPICseq_all_merged_bams_output/coverage.2k.matrix.merged.by.60.txt", sep = "\t")
Cov_2k_gr60 <- Cov_2k_gr60[, grep("NDR_all", colnames(Cov_2k_gr60))]
Cov_2k_vector_gr60<-as.vector(unlist(Cov_2k_gr60))

Cov_2k_gr80 <- read.csv("~/rds/hpc-work/EPICseq_all_merged_bams_output/coverage.2k.matrix.merged.by.80.txt", sep = "\t")
Cov_2k_gr80 <- Cov_2k_gr80[, grep("NDR_all", colnames(Cov_2k_gr80))]
Cov_2k_vector_gr80<-as.vector(unlist(Cov_2k_gr80))

Cov_2k_gr100 <- read.csv("~/rds/hpc-work/EPICseq_all_merged_bams_output/coverage.2k.matrix.merged.by.100.txt", sep = "\t")
Cov_2k_gr100 <- Cov_2k_gr100[, grep("NDR_all", colnames(Cov_2k_gr100))]
Cov_2k_vector_gr100<-as.vector(unlist(Cov_2k_gr100))


Cov_2k_list <- list(Group_by_0 = Cov_2k_vector_gr0,
                    Group_by_20 = Cov_2k_vector_gr20,
                    Group_by_40 = Cov_2k_vector_gr40,
                    Group_by_60 = Cov_2k_vector_gr60,
                    Group_by_80 = Cov_2k_vector_gr80,
                    Group_by_100 = Cov_2k_vector_gr100)


pdf("~/rds/hpc-work/output_Rplots/mergedvisits_coverage_boxplot.pdf", width = 7, height = 5)
boxplot(Cov_2k_list, names = seq(0,100,20), col = "lightblue", main = "Coverage over 2kbp region surrounding each TSS (merged samples)",
        xlab = "Number of genes per group",
        ylab = "Number of fragments in 2kbp region",
        ylim = c(0, 5000))
dev.off()



# ****** Investigating NDR coverage : per visit files ****** #
NDR_cov <- read.csv("~/rds/hpc-work/Output_Main/coverage.ndr.matrix.merged.by.100.txt", sep = "\t")
NDR_cov <- NDR_cov[, grep("NDR_all", colnames(NDR_cov))]

NDR_cov_vector<-as.vector(unlist(NDR_cov))
median(NDR_cov_vector) 
range(NDR_cov_vector) 

# ****** Investigating 2k coverage : per visit files ****** #

Cov_2k_gr0 <- read.csv("~/rds/hpc-work/Output_Main/coverage.2k.matrix.merged.by.0.txt", sep = "\t")
Cov_2k_gr0 <- Cov_2k_gr0[, grep("NDR_all", colnames(Cov_2k_gr0))]
Cov_2k_vector_gr0<-as.vector(unlist(Cov_2k_gr0))

Cov_2k_gr20 <- read.csv("~/rds/hpc-work/Output_Main/coverage.2k.matrix.merged.by.20.txt", sep = "\t")
Cov_2k_gr20 <- Cov_2k_gr20[, grep("NDR_all", colnames(Cov_2k_gr20))]
Cov_2k_vector_gr20<-as.vector(unlist(Cov_2k_gr20))

Cov_2k_gr40 <- read.csv("~/rds/hpc-work/Output_Main/coverage.2k.matrix.merged.by.40.txt", sep = "\t")
Cov_2k_gr40 <- Cov_2k_gr40[, grep("NDR_all", colnames(Cov_2k_gr40))]
Cov_2k_vector_gr40<-as.vector(unlist(Cov_2k_gr40))

Cov_2k_gr60 <- read.csv("~/rds/hpc-work/Output_Main/coverage.2k.matrix.merged.by.60.txt", sep = "\t")
Cov_2k_gr60 <- Cov_2k_gr60[, grep("NDR_all", colnames(Cov_2k_gr60))]
Cov_2k_vector_gr60<-as.vector(unlist(Cov_2k_gr60))

Cov_2k_gr80 <- read.csv("~/rds/hpc-work/Output_Main/coverage.2k.matrix.merged.by.80.txt", sep = "\t")
Cov_2k_gr80 <- Cov_2k_gr80[, grep("NDR_all", colnames(Cov_2k_gr80))]
Cov_2k_vector_gr80<-as.vector(unlist(Cov_2k_gr80))

Cov_2k_gr100 <- read.csv("~/rds/hpc-work/Output_Main/coverage.2k.matrix.merged.by.100.txt", sep = "\t")
Cov_2k_gr100 <- Cov_2k_gr100[, grep("NDR_all", colnames(Cov_2k_gr100))]
Cov_2k_vector_gr100<-as.vector(unlist(Cov_2k_gr100))


Cov_2k_list <- list(Group_by_0 = Cov_2k_vector_gr0,
                    Group_by_20 = Cov_2k_vector_gr20,
                    Group_by_40 = Cov_2k_vector_gr40,
                    Group_by_60 = Cov_2k_vector_gr60,
                    Group_by_80 = Cov_2k_vector_gr80,
                    Group_by_100 = Cov_2k_vector_gr100)


pdf("~/rds/hpc-work/output_Rplots/pervisit_coverage_boxplot.pdf", width = 7, height = 5)
boxplot(Cov_2k_list, names = seq(0,100,20), col = "lightblue", main = "Coverage over 2kbp region surrounding each TSS (per visit samples)",
        xlab = "Number of genes per group",
        ylab = "Number of fragments in 2kbp region",
        ylim = c(0,5000))
dev.off()


patient_subtype_data <- read.csv("~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/DETECT_PublicID_ER.Her2.txt", sep = "\t")


# ********* MERGED SAMPLES : Inferred GEP grouped by 20 ********* #

GEP_data <- read.csv("~/rds/hpc-work/GEP_output_all_merged_bams/Grouped_by_20/EPICSeqInferredExpressionValues.txt", sep = "\t")
((length(which(GEP_data$inferredGEP == "-Inf")) + length(which(is.na(GEP_data$inferredGEP)))) / nrow(GEP_data)) * 100
# 20 = 44%

pdf("~/rds/hpc-work/output_Rplots/Dist_of_inferredGEP_plots/All_merged_bams/GEP_dist_plot_merged_gr20.pdf", width = 7, height = 5)
hist(GEP_data$inferredGEP, 
     xlab = "Inferred gene expression", 
     ylab = "Gene group count", 
     col = "deepskyblue3",
     breaks = 30,
     main = "Inferred GEP distribution for merged files (grouped per 20 genes)")
box(col = "black", lwd = 1)  # 'col' sets the border color, 'lwd' sets the border line width
dev.off()

# ***** ESR1 gene ***** #
ESR1_data <- GEP_data[grepl("ESR1", GEP_data$gene),]

ESR1_data$patient <- NA
ESR1_data$Subtype <- NA
for (row in 1:nrow(ESR1_data)){
  patient <- gsub("merged_", "", ESR1_data$Sample[row])
  patient <- gsub(".sorted", "", patient)
  ESR1_data$patient[row] <- patient
  subtype <- patient_subtype_data[which(patient_subtype_data$Public.ID == patient), "ER.status"]
  ESR1_data$Subtype[row] <- unique(subtype)
}

plot5 <- ggplot(ESR1_data, aes(x = Subtype, y = inferredGEP, fill = Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Inferred GEP") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold", size = 14), 
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
        panel.border = element_rect(color = "black", size = 1)) +
  ggtitle("No significant difference in inferred ESR1 expression between ER+ and ER- subtypes") +
  scale_y_continuous(limits = c(-5, 0))

ggsave(filename = paste("~/rds/hpc-work/output_Rplots/Dist_of_inferredGEP_plots/All_merged_bams/ESR1_mergedby20.png", sep = ""),
       plot5, width = 9, height = 5, device = png)


kruskal.test(formula = inferredGEP ~ Subtype, data = ESR1_data) # Kruskal-Wallis chi-squared = 0.67818, df = 1, p-value = 0.4102

# ***** ERBB2 gene ***** #
ERBB2_data <- GEP_data[grepl("ERBB2", GEP_data$gene),]

ERBB2_data$patient <- NA
ERBB2_data$Subtype <- NA
for (row in 1:nrow(ERBB2_data)){
  patient <- gsub("merged_", "", ERBB2_data$Sample[row])
  patient <- gsub(".sorted", "", patient)
  ERBB2_data$patient[row] <- patient
  subtype <- patient_subtype_data[which(patient_subtype_data$Public.ID == patient), "Her2.status"]
  ERBB2_data$Subtype[row] <- unique(subtype)
}

plot6 <- ggplot(ERBB2_data, aes(x = Subtype, y = inferredGEP, fill = Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Inferred GEP") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold", size = 14), 
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
        panel.border = element_rect(color = "black", size = 1)) +
  ggtitle("No significant difference in inferred ERBB2 expression between HER2+ and HER2- subtypes") +
  scale_y_continuous(limits = c(-4, 2))

ggsave(filename = paste("~/rds/hpc-work/output_Rplots/Dist_of_inferredGEP_plots/All_merged_bams/ERBB2_mergedby20.png", sep = ""),
       plot6, width = 9, height = 5, device = png)


kruskal.test(formula = inferredGEP ~ Subtype, data = ERBB2_data) # Kruskal-Wallis chi-squared = 0.043569, df = 1, p-value = 0.8347


# ********* PER VISIT SAMPLES : Inferred GEP grouped by 20 ********* #

GEP_data <- read.csv("~/rds/hpc-work/GEP_output_all_BAMS/mergedby20/EPICSeqInferredExpressionValues.txt", sep = "\t")
(length(which(GEP_data$inferredGEP == "-Inf")) + length(which(is.na(GEP_data$inferredGEP)))) / nrow(GEP_data)

# 0 = 98%, 20 = 98%, 40 = 98%

pdf("~/rds/hpc-work/output_Rplots/Dist_of_inferredGEP_plots/Per_visit_bams/GEP_dist_plot_pervisit_gr20.pdf", width = 7, height = 5)
hist(GEP_data$inferredGEP, 
     xlab = "Inferred gene expression", 
     ylab = "Gene group count", 
     col = "deepskyblue3",
     breaks = 30,
     main = "Inferred GEP distribution for per visit files (grouped per 20 genes)")
box(col = "black", lwd = 1)  # 'col' sets the border color, 'lwd' sets the border line width
dev.off()

# ***** ESR1 gene ***** #
ESR1_data <- GEP_data[grepl("ESR1", GEP_data$gene),]

ESR1_data$patient <- NA
ESR1_data$Subtype <- NA
for (row in 1:nrow(ESR1_data)){
  patient <- gsub("merged_", "", ESR1_data$Sample[row])
  patient <- gsub(".sorted", "", patient)
  ESR1_data$patient[row] <- patient
  subtype <- patient_subtype_data[which(patient_subtype_data$Public.ID == patient), "ER.status"]
  ESR1_data$Subtype[row] <- unique(subtype)
}

plot5 <- ggplot(ESR1_data, aes(x = Subtype, y = inferredGEP, fill = Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Inferred GEP") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"), 
        axis.title = element_text(size = 12)) +
  ggtitle("Barplot showing difference in inferred ESR1 gene expression between subtypes")

ggsave(filename = paste("~/rds/hpc-work/output_Rplots/Dist_of_inferredGEP_plots/All_merged_bams/ESR1_mergedby20.png", sep = ""),
       plot5, width = 9, height = 5, device = png)


# ***** ERBB2 gene ***** #
ERBB2_data <- GEP_data[grepl("ERBB2", GEP_data$gene),]

ERBB2_data$patient <- NA
ERBB2_data$Subtype <- NA
for (row in 1:nrow(ERBB2_data)){
  patient <- gsub("merged_", "", ERBB2_data$Sample[row])
  patient <- gsub(".sorted", "", patient)
  ERBB2_data$patient[row] <- patient
  subtype <- patient_subtype_data[which(patient_subtype_data$Public.ID == patient), "Her2.status"]
  ERBB2_data$Subtype[row] <- unique(subtype)
}

plot6 <- ggplot(ERBB2_data, aes(x = Subtype, y = inferredGEP, fill = Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Inferred GEP") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,vjust = 3, face = "bold"), 
        axis.title = element_text(size = 12)) +
  ggtitle("Barplot showing difference in inferred ERBB2 gene expression between subtypes")

ggsave(filename = paste("~/rds/hpc-work/output_Rplots/Dist_of_inferredGEP_plots/All_merged_bams/ERBB2_mergedby20.png", sep = ""),
       plot6, width = 9, height = 5, device = png)


