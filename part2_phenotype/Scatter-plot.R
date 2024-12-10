###相关性###
# 选择要比较的列名
column_name <- "TotalIPA_2015"  # 替换为实际的列名


# 原始数据的散点图
plot(pheno_matrix[, column_name], pheno_n[, column_name],
     xlab = "Original Data", ylab = "Standardized Data",
     main = paste("Scatter Plot of", column_name))


# 添加回归线
abline(lm(pheno_n[, column_name] ~ pheno_matrix[, column_name]), col = "red")
