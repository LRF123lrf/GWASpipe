###协变量文件的准备格式：FID  IID  截距1  PC1 PC2 PC3...PCn###


### 读取文件 ###
tfam <- read.table("emmax_in.tfam", header = FALSE)
pca <- read.table("pca_output.eigenvec", header = FALSE)

### 提取前两列和 PCA 的第3到11列 ###
tfam_cols <- tfam[, 1:2]
pca_cols <- pca[, 3:11]

### 创建全为1的第三列 ###
third_col <- rep(1, nrow(tfam))

### 合并所有列 ###
final_data <- cbind(tfam_cols, third_col, pca_cols)

### 保存为 emmax.cov.txt 文件 ###
write.table(final_data, file = "emmax.cov.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
