# library
library(BiocGenerics)

# set info
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
    
    # filter
    gtf_temp <- gtf_temp[gtf_temp@seqnames == chr]
    gtf_temp <- gtf_temp[order(start(gtf_temp))]
    # cat(length(gtf_temp), "\n")

    # make data frame
    if (i == 2) {
        gtf_df <- Repitools::annoGR2DF(gtf_temp)
    } else {
        gtf_df$cov <- gtf_temp$cov
        gtf_df$FPKM <- gtf_temp$FPKM
        gtf_df$TPM <- gtf_temp$TPM
    }
    colnames(gtf_df)[which(colnames(gtf_df)=="cov")] <- paste("cov", "_", gtf_name)
    gtf_df$cov <- NULL
    colnames(gtf_df)[which(colnames(gtf_df)=="FPKM")] <- paste("FPKM", "_", gtf_name)
    gtf_df$FPKM <- NULL
    colnames(gtf_df)[which(colnames(gtf_df)=="TPM")] <- paste("TPM", "_", gtf_name)
    gtf_df$TPM <- NULL
    
}
cat("--------------loop end--------------", "\n")

# save csv
cat("-------------writing csv-------------", "\n")
write.csv(gtf_df, file.path(path, paste0(output_name, ".csv")))