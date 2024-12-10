###绘图###

#绘制处理前的表型正态分布图
data <- read.csv("pheno.csv", header = TRUE, row.names = 1)
# 将宽格式数据转换为长格式
data_long <- pivot_longer(data, cols = everything(), names_to = "Phenotype", values_to = "Value")
# 绘制正态分布图
ggplot(data_long, aes(x = Value, fill = Phenotype)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Phenotype, scales = "free") +  # 每个表型单独绘制
  labs(title = "Normal Distribution of Phenotypes", x = "Value", y = "Density") +
  theme_minimal()

#绘制处理后的表型正态分布图
data <- read.csv("pheno_standardized.csv", header = TRUE, row.names = 1)
# 将宽格式数据转换为长格式
data_long <- pivot_longer(data, cols = everything(), names_to = "Phenotype", values_to = "Value")
# 绘制正态分布图
ggplot(data_long, aes(x = Value, fill = Phenotype)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Phenotype, scales = "free") +  # 每个表型单独绘制
  labs(title = "Normal Distribution of Phenotypes", x = "Value", y = "Density") +
  theme_minimal()

