
# Read in old files
Ch37_genes <- read.csv("~/rds/hpc-work/Code_and_prior_data/priordata/all.tss.genes.canonical.ensembl75.old.txt", sep = "\t")

Ch38_genes <- read.csv("BioMart_GeneInfo3.txt") 
Ch38_genes <- Ch38_genes[!duplicated(Ch38_genes$Gene.name), ]

Reference_old <- read.csv("~/rds/hpc-work/Code_and_prior_data/priordata/reference.pbmc.RNA.TPM.old.txt", sep = "\t")
Reference_exome_old <- read.csv("~/rds/hpc-work/Code_and_prior_data/priordata/reference.pbmc.RNA.TPM.ExomeSpace.old.txt", sep = "\t")


Ch38_genes$TSS <- NA
for (row in 1:nrow(Ch38_genes)){
  if(Ch38_genes[row, "Strand"] == 1){
    Ch38_genes$TSS[row] <- Ch38_genes[row, "Gene.start..bp."]
  } else {
    Ch38_genes$TSS[row] <- Ch38_genes[row, "Gene.end..bp."]
  }
}

Ch37_genes$TSS_38 <- NA
Ch37_genes$gc <- NA
for (row in 1:nrow(Ch37_genes)) {
  Gene_name <- Ch37_genes$Gene.Symbol[row]
  matching_indices <- which(Ch38_genes$Gene.name == Gene_name)
  matching_indices_alt <- which(Ch38_genes$Gene.Synonym == Gene_name)
  
  if (length(matching_indices) > 0) {
    TSS_38 <- unique(Ch38_genes$TSS[matching_indices])
    GC <- unique(Ch38_genes$Gene...GC.content[matching_indices])
    Ch37_genes$TSS_38[row] <- TSS_38
    Ch37_genes$gc[row] <- GC/100
  } else if (length(matching_indices) == 0 && length(matching_indices_alt) > 0){
    TSS_38 <- unique(Ch38_genes$TSS[matching_indices_alt])
    GC <- unique(Ch38_genes$Gene...GC.content[matching_indices_alt])
    Ch37_genes$TSS_38[row] <- TSS_38
    Ch37_genes$gc[row] <- GC/100
  } else {
    Ch37_genes$TSS_38[row] <- NA
    Ch37_genes$gc[row] <- NA
  }
}

# Update reference genome file
Reference_old$Ch38_gene_start <- NA
Reference_old$Ch38_gene_end <- NA

for (row in 1:nrow(Reference_old)) {
  Gene_name <- Reference_old$geneName[row]
  matching_indices <- which(Ch38_genes$Gene.name == Gene_name)
  matching_indices_alt <- which(Ch38_genes$Gene.Synonym == Gene_name)
  
  if (length(matching_indices) > 0) {
    New_start <- unique(Ch38_genes$Gene.start..bp.[matching_indices])
    New_end <- unique(Ch38_genes$Gene.end..bp.[matching_indices])
    Reference_old$Ch38_gene_start[row] <- New_start
    Reference_old$Ch38_gene_end[row] <- New_end
  } else if (length(matching_indices) == 0 && length(matching_indices_alt) > 0){
    New_start <- unique(Ch38_genes$Gene.start..bp.[matching_indices_alt])
    New_end <- unique(Ch38_genes$Gene.end..bp.[matching_indices_alt])
    Reference_old$Ch38_gene_start[row] <- New_start
    Reference_old$Ch38_gene_end[row] <- New_end
  } else {
    Reference_old$Ch38_gene_start[row] <- NA
    Reference_old$Ch38_gene_end[row] <- NA
  }
}

# Update gene start and end column to produce updated reference file
Reference_updated <- Reference_old
Reference_updated$start <- Reference_updated$Ch38_gene_start
Reference_updated$end <- Reference_updated$Ch38_gene_end
Reference_updated$Ch38_gene_start <- NULL
Reference_updated$Ch38_gene_end <- NULL

Reference_NAs <- which(is.na(Reference_updated$start) == TRUE)
Reference_updated <- Reference_updated[-c(Reference_NAs),]


# Update bed file
bed_file <- read.csv("~/rds/hpc-work/Code_and_prior_data/priordata/ensembl75.wgEncodeDukeMapabilityUniqueness35bp.old.bed", sep = "\t", header = F)
colnames(bed_file) <- c("chromosome", "minus_1kbp", "plus_1kbp", "ID", "uniqueness")

bed_file_updated <- bed_file
bed_file_updated$TSS_37 <- (bed_file$minus_1kbp + bed_file$plus_1kbp)/2


# Use their BED file and find their TSS by taking middle of their genomic coordinates
# Use the Ch37_genes file to find our corresponding TSS
# Update BED file with -1kb and +1kb from our TSS
for (row in 1:nrow(bed_file_updated)) {
  TSS_37 <- bed_file_updated$TSS_37[row]
  TSS_38 <- Ch37_genes[which(Ch37_genes$TSS == TSS_37), "TSS_38"]
  
  if (length(TSS_38) == 1 && !is.na(TSS_38)) {
    bed_file_updated$minus_1kbp[row] <- TSS_38 - 1000
    bed_file_updated$plus_1kbp[row] <- TSS_38 + 1000
  } else if (length(TSS_38) == 1 && is.na(TSS_38)) {
    bed_file_updated$minus_1kbp[row] <- "TSS_38_is_NA"
    bed_file_updated$plus_1kbp[row] <- "TSS_38_is_NA"
  } else if (length(TSS_38) > 1 && !is.na(TSS_38)) {
    bed_file_updated$minus_1kbp[row] <- "ambiguous"
    bed_file_updated$plus_1kbp[row] <- "ambiguous"
  } else if (length(TSS_38) > 1 && is.na(TSS_38)) {
    bed_file_updated$minus_1kbp[row] <- "ambiguous_and_TSS_38_NA"
    bed_file_updated$plus_1kbp[row] <- "ambiguous_and_TSS_38_NA"
  }
}

# Remove ambiguous and NAs from BED file
TSS_38_is_NA <- which(bed_file_updated$minus_1kbp == "TSS_38_is_NA")
Ambiguous <- which(bed_file_updated$minus_1kbp == "ambiguous")
Ambiguous_and_TSS_38_NA <- which(bed_file_updated$minus_1kbp == "ambiguous_and_TSS_38_NA")
bed_file_updated <- bed_file_updated[-c(TSS_38_is_NA, Ambiguous, Ambiguous_and_TSS_38_NA),]

# Update ID column in BED file
for (row in 1:nrow(bed_file_updated)){
  chromosome <- bed_file_updated$chromosome[row]
  minus_1kb <- bed_file_updated$minus_1kbp[row]
  plus_1kb <- bed_file_updated$plus_1kbp[row]
  bed_file_updated$ID[row] <- paste(chromosome, "_", minus_1kb, "_", plus_1kb, sep = "")
}

bed_file_updated$TSS_37 <- NULL

# Remove NAs from TSS info file
TSS_38_NAs <- which(is.na(Ch37_genes$TSS_38) == TRUE)
Ch37_genes <- Ch37_genes[-c(TSS_38_NAs),]
# Update TSS col
Ch37_genes$TSS <- Ch37_genes$TSS_38
Ch37_genes$TSS_38 <- NULL

# Create GC dataframe
GC_df <- data.frame("gene" = Ch37_genes$Gene.Symbol, "gc" = Ch37_genes$gc)
GC_df <- GC_df[order(GC_df$gene), ]

Ch37_genes$gc<-NULL
colnames(Ch37_genes)[3] <- "Gene-Symbol"


# Write files
write.table(Ch37_genes, file = "~/rds/hpc-work/Code_and_prior_data/priordata/all.tss.genes.canonical.ensembl75.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bed_file_updated, file = "~/rds/hpc-work/Code_and_prior_data/priordata/ensembl75.wgEncodeDukeMapabilityUniqueness35bp.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(GC_df, file = "~/rds/hpc-work/Code_and_prior_data/priordata/ensembl75.TSS_GC_aggregate.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Reference_updated, file = "~/rds/hpc-work/Code_and_prior_data/priordata/reference.pbmc.RNA.TPM.txt", sep = "\t", row.names = FALSE, quote = FALSE)






