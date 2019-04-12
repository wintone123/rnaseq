# package
library(dplyr)
library(tidyr)
library(parallel)

# functions
AlterRatio <- function(name) {
    # arrange
    data <- filter(input, ref_gene_name == name)
    data <- gather(data, c(3:ncol(data)), key = "sample", value = "value")
    data <- data[complete.cases(data$value),]
    data <- separate(data, "sample", c("type", "sex", "pop", "num"), sep = "_")
    data <- unite(data, "name", c("sex", "pop"), sep = "_")
    
    # calculate
    data <- group_by(data, ref_gene_name, transcript_id, type, name)
    data <- summarise(data, mean = mean(value), sd = sd(value))
    data$ratio <- NA
    for (name in unique(data$name)) {
        for (type in unique(data$type)) {
            sub_data <- data[data$name == name & data$type == type,]
            data[data$name == name & data$type == type,]$ratio <- sub_data$mean / sum(sub_data$mean)
        }
    }   
      
    return(data)
}

# settings
path = "/mnt/d/RNASeq/samples/chrX_data/chrX_data/HISAT/ballgown"
imput_csv = "output.csv"
output_csv = "altersplice.csv"

# load file
cat("------------loading file-------------", "\n")
input = read.delim(file.path(path, imput_csv), sep = ",")

# arrange
input <- input[complete.cases(input$ref_gene_name),]
input <- filter(input, type == "transcript")
input <- select(input, c("ref_gene_name", "transcript_id", starts_with("FPKM"), starts_with("TPM")))

# 2 thread processing
cat("-------------processing--------------", "\n")
cl <- makeForkCluster(2)
clusterExport(cl, "input")
temp <- parLapply(cl, unique(input$ref_gene_name), AlterRatio)
output <- Reduce("rbind", temp)
stopCluster(cl)

# arrange
output <- arrange(output, ref_gene_name, name, type, transcript_id)
output <- output[,c(1:5,7,6)]

# output
cat("------------writing csv--------------", "\n")
write.table(output, file.path(path, output_csv), sep = ",", row.names = FALSE)