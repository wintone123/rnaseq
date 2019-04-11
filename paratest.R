# package
library(parallel)
library(tidyr)
library(dplyr)

# functions
fun <- function(i) {
    k <- 1
    for (j in data$a[i]:data$b[i]) {
        k <- k * j
    }
    temp <- data[i,]
    temp$out[1] <- k
    return(temp)
}

# settings
path = "/mnt/d/Desktop"
imput_csv = "input.csv"
output_csv = "output.csv"

# load file
# cat("------------loading file-------------", "\n")
data = read.delim(file.path(path, imput_csv), sep = ",")

# filtering

# make new columns
data$out <- rep(NA, nrow(data))

# processing
cat("-------------processing--------------", "\n")
list <- c(1:nrow(data))
print(Sys.time())
cl <- makeForkCluster(2)
clusterExport(cl, "data")
temp <- parLapply(cl, list, fun)
output <- do.call("rbind", temp)
stopCluster(cl)
# for (i in list) {
#     k <- 1
#     for (j in data$a[i]:data$b[i]) {
#         k <- k * j
#     }
#     data$out[i] <- k
# }
print(Sys.time())

# output
cat("------------writing csv--------------", "\n")
write.table(data, file.path(path, output_csv), sep = ",", row.names = FALSE)