# 读取 CSV 文件
pheno_matrix <- read.csv("pheno_log.csv", header = TRUE, row.names = 1)
# 如果需要将数据框转换为矩阵
pheno_matrix <- as.matrix(pheno_matrix)


# 假设 rankTransPheno 函数已经定义好
rankTransPheno <- function(pheno, para_c) {
  if (length(pheno) == 0) {
    return(pheno)
  }
  # 计算排名并映射到正态分布
  pheno <- qnorm((rank(pheno, na.last = "keep") - para_c) / (length(pheno) - 2 * para_c + 1))
  return(pheno)
}

# 定义 rankTransPheno_multi 函数
rankTransPheno_multi <- function(pheno_matrix, para_c) {
  pheno_norm <- apply(pheno_matrix, 2, function(phei) {
    # 获取列名
    namei <- names(phei)
    # 提取缺失值和非缺失值
    phe.na <- phei[is.na(phei)]
    phei <- phei[!is.na(phei)]
    # 进行排名转换
    phei <- rankTransPheno(phei, para_c)
    # 恢复缺失值的位置
    phei <- c(phei, rep(NA, length(phe.na)))[namei]
    return(phei)
  })
  return(pheno_norm)
}

# 应用函数
pheno_n <- rankTransPheno_multi(pheno_matrix, 0.5)

# 查看结果
print(pheno_n)
# 保存标准化后的数据为 CSV 文件
write.csv(pheno_n, file = "pheno_standardized.csv", row.names = TRUE)
