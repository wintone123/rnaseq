# library
library(BiocGenerics)
library(dplyr)

# settings
path = "/mnt/d/RNASeq/samples/chrX_data/chrX_data/HISAT/ballgown"
output_name = "output"
chr = "chrX"

# load path
cat("------------loading path------------", "\n")
dic_list <- list.dirs(path)

# process
cat("-------------loop start-------------", "\n")
for(i in 2:length(dic_list)) {
    dic <- dic_list[i]
    gtf_name <- tail(strsplit(dic, split = "/", fixed = TRUE)[[1]], n = 1)
    cat("-------processing",gtf_name,"--------", "\n")

    # load gft
    gtf_temp <- rtracklayer::import.gff2(file.path(path, gtf_name, paste0(gtf_name, ".gtf")))
    
    # filtering
    gtf_temp <- gtf_temp[gtf_temp@seqnames == chr]
    gtf_temp <- gtf_temp[order(start(gtf_temp))]
    # cat(length(gtf_temp), "\n")

    # make data frame
    if (i == 2) {
        gtf_df <- Repitools::annoGR2DF(gtf_temp)
        gtf_df <- gtf_df[,c(1,14:16,2:13)]
        gtf_df$cov <- NULL
    } else {
        # gtf_df$cov <- gtf_temp$cov
        gtf_df$FPKM <- gtf_temp$FPKM
        gtf_df$TPM <- gtf_temp$TPM
    }
    # colnames(gtf_df)[which(colnames(gtf_df)=="cov")] <- paste0("cov", "_", gtf_name)
    # gtf_df$cov <- NULL
    colnames(gtf_df)[which(colnames(gtf_df)=="FPKM")] <- paste0("FPKM", "_", gtf_name)
    gtf_df$FPKM <- NULL
    colnames(gtf_df)[which(colnames(gtf_df)=="TPM")] <- paste0("TPM", "_", gtf_name)
    gtf_df$TPM <- NULL
}
cat("--------------loop end--------------", "\n")

# arranging
gtf_df <- filter(gtf_df, source == "StringTie")
gtf_df$source <- NULL
gtf_df$phase <- NULL
gtf_df$type <- NULL
gtf_df$gene_name <- NULL

# save csv
cat("-------------writing csv-------------", "\n")
write.csv(gtf_df, file.path(path, paste0(output_name, ".csv")))