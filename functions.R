# alternative splices map
AlterMap <- function(input) {
    alter <- GRanges(seqnames = input$chr,
                     ranges = IRanges(start = input$start, end = input$end),
                     strand = Rle("*"),
                     ref_gene_name = input$ref_gene_name,
                     type = input$type,
                     gene = input$gene_id,
                     transcript = input$transcript_id)
    
    # filter 
    alter <- alter[alter$type == "exon"]
    alter <- alter[order(alter$transcript)]
        
    # plot
    alter_track <- GeneRegionTrack(alter, name = as.character(alter$ref_gene_name[1]),
                                   transcriptAnnotation = "transcript")
    plotTracks(alter_track)
}

# alternative splices ratio
AlterRatio <- function(input) {
    # filter
    input <- filter(input, type == "transcript") 
    data <- select(input, c("ref_gene_name", "transcript_id", starts_with("FPKM"), starts_with("TPM")))
    data <- gather(data, c(3:ncol(data)), key = "sample", value = "value")
    data <- separate(data, "sample", c("type", "name", "num"), sep = "_")
    
    # calculate
    data <- group_by(data, ref_gene_name, transcript_id, type, name)
    data <- summarise(data, mean = mean(value))
    data$ratio <- NA
    for (name in unique(data$name)) {
        sub_data <- data[data$name == name,]
        data[data$name == name,]$ratio <- sub_data$mean / ifelse(sub_data$type == "FPKM", sum(sub_data[sub_data$type == "FPKM",]$mean), sum(sub_data[sub_data$type == "TPM",]$mean))
    }
    data <- arrange(data, ref_gene_name, name, type, transcript_id)
    
    return(data)
}

# 