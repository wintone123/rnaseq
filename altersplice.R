# package
library(dplyr)
library(tidyr)

# functions
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
        for (type in unique(data[data$name == name,]$type)) {
            sub_data <- data[data$name == name & data$type == type,]
            data[data$name == name & data$type == type,]$ratio <- sub_data$mean / sum(sub_data$mean)
        }
    }
      
    return(data)
}

# settings
path = "/mnt/d/RNASeq/samples/chrX_data/chrX_data/HISAT/ballgown"
imput_csv = "output.csv"
output_name = "altersplice"

# load file
cat("------------loading file-------------", "\n")
input = read.delim(file.path(path, imput_csv), sep = ",")

# filtering
input <- input[complete.cases(input$ref_gene_name),]

# processing
cat("-------------processing--------------", "\n")
name_list <- unique(input$ref_gene_name)
j <- 0
k <- 0
for (i in 1:length(name_list)) {
    if (i == 1) {
        output <- AlterRatio(filter(input, ref_gene_name == name_list[i]))
    } else {
        output <- rbind(output, AlterRatio(filter(input, ref_gene_name == name_list[i])))
    }

    # progress bar
    if ((j / length(name_list) * 100) %/% 1 == k) { 
        cat("[", paste(rep("#", (31*k/100) %/% 1),collapse = ""), 
        paste(rep("_", 32-(32*k/100) %/% 1),collapse = ""), k, "%]","\r", sep = "")
        k <- k + 1
    }
    j <- j + 1
}
output <- arrange(output, ref_gene_name, name, type, transcript_id)

# output
cat("-------------writing csv-------------", "\n")
write.csv(output, file.path(path, paste0(output_name, ".csv")))