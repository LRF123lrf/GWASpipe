# 指定表型目录和输出目录
pheno_dir <- "./emmax_pheno/"
output_dir <- "./emmax_scripts/"

# 获取表型文件列表
pheno_files <- list.files(pheno_dir, pattern = "_Pheno.txt")

# 遍历每个表型文件并生成相应的脚本
for (pheno_file in pheno_files) {
  # 提取表型名称（去掉_Pheno.txt部分）
  pheno_name <- gsub("_Pheno.txt", "", pheno_file)
  
  # 定义脚本内容
  script_content <- paste0(
    "#BSUB -J ", pheno_name, "_kq.sh\n",
    "#BSUB -n 5\n",
    "#BSUB -R \"span[hosts=1]\"\n",
    "#BSUB -o ", pheno_name, "_kq.sh.%J.out\n",
    "#BSUB -e ", pheno_name, "_kq.sh.%J.err\n",
    "#BSUB -q normal\n\n",

    "/public/home/rfluo/01.software/emmax-intel64 -v -d 10 -t emmax_in -o output/", pheno_name, ".kq -p ", pheno_dir, pheno_file, " -k emmax_in.BN.kinf -c emmax.cov.txt\n"
  )
  
  # 将脚本保存为文件
  script_file <- paste0(output_dir, pheno_name, "_kq.sh")
  writeLines(script_content, script_file)
  
  # 打印生成的脚本路径以确认
  print(paste("Generated script:", script_file))
}
