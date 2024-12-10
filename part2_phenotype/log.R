# 读取CSV文件
data <- read.csv("pheno.csv")
# 查看数据的结构
str(data)
data_log2 <- data
data_log2[sapply(data, is.numeric)] <- log2(data[sapply(data, is.numeric)] + 1)
head(data_log2)
write.csv(data_log2, "pheno_log.csv", row.names = FALSE)
