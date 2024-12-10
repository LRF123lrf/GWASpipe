###表型文件的处理格式为三列：FID  IID  pheno-data###

### 读取 .tfam 文件 ###
tfam <- read.table("emmax_in.tfam", header = FALSE)
# 读取 pheno_n_standardized.csv 文件
pheno <- read.csv("pheno_standardized.csv", header = TRUE)

### 创建输出目录 ###
output_dir <- "./emmax_pheno"

### 提取 tfam 文件的前两列 ###
tfam_cols <- tfam[, 1:2]
### 遍历 pheno 的每一列（除了第一列）###
for (i in 2:ncol(pheno)) {
  ### 提取 pheno 的当前列 ###
  pheno_col <- pheno[, i]
  
  ### 合并 tfam 的前两列和 pheno 的当前列 ###
  combined_data <- cbind(tfam_cols, pheno_col)
  
  ### 获取当前列的列名，并创建输出文件名 ###
  pheno_name <- colnames(pheno)[i]
  output_filename <- paste0(output_dir, "/", pheno_name, "_Pheno.txt")
  
  ### 输出合并后的数据到文件 ###
  write.table(combined_data, file = output_filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # 打印当前生成的文件名，确认输出
  print(paste("Generated file:", output_filename))
}
