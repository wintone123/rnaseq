# output alternative splices
alternatives <- function(gene_name, input) {
    input <- input[complete.cases(input$ref_gene_name),]
    alter <- GRanges(seqnames = input$chr,
                     ranges = IRanges(start = input$start, end = input$end),
                     strand = Rle("*"),
                     ref_gene_name = input$ref_gene_name,
                     type = input$type,
                     gene = input$gene_id,
                     transcript = input$transcript_id)
    alter <- alter[alter$ref_gene_name == gene_name]
    alter <- alter[order(alter$transcript)]
    alter <- alter[alter$type == "exon"]
    alter_track <- GeneRegionTrack(alter, name = gene_name,
                                   transcriptAnnotation = "transcript")
    plotTracks(alter_track)
}