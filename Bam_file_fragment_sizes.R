if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("csaw")

library("csaw")

directory <- "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/ctDNA/data/DETECT/bam_files/Paired_End"

bam_files <- list.files(directory, pattern = "\\.bam$", full.names = TRUE)

# For each BAM file, record sizes of each fragment and write to a text file
for (bam_file in bam_files) {
  sizes <- getPESizes(bam_file, param=readParam(pe="both"))$sizes
  file_name<-basename(bam_file)
  file_name <- sub("\\.bam$", "", file_name)
  sample <- sub("^P", "", file_name)
  write.table(sizes, paste("~/rds/hpc-work/Read_length_outputfiles/", sample, ".txt",  sep = ""), row.names = FALSE, col.names = FALSE)
}




