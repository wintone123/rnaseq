# package
library(dplyr)
library(tidyr)

# functions

# settings
path = "/mnt/d/RNASeq/samples/chrX_data/chrX_data/HISAT/ballgown"
imput_csv = "output.csv"
output_csv = "total.csv"

# load file
cat("------------loading file-------------", "\n")
input = read.delim(file.path(path, imput_csv), sep = ",")

# arrange
input <- input[complete.cases(input$ref_gene_name),]
input <- filter(input, type == "transcript")

# processing
cat("-------------processing--------------", "\n")
output <- gather(input, c(11:34), key = "sample", value = "value")
output <- output[complete.cases(output$value),]
output <- group_by(output, ref_gene_name, gene_id, sample)
output <- summarise(output, sum = sum(value))
output <- separate(output, sample, c("type", "sex", "pop", "num"))
output <- unite(output, "name", c("sex", "pop"), sep = "_")
output <- group_by(output, ref_gene_name, type, name)
output <- summarise(output, mean = mean(sum), sd = sd(sum))

# output
cat("------------writing csv--------------", "\n")
write.table(output, file.path(path, output_csv), sep = ",", row.names = FALSE)