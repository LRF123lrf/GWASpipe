
rule all:
    input:
        "result/04.anno/SNP.vcf",
        "result/05.plink/ld.prune.in.bed",
        "result/05.plink/ld.prune.in.bim",
        "result/05.plink/ld.prune.in.fam"

rule filter_maf_gene:
    input:
        "result/04.anno/anno.vcf.gz"
    output:
        "result/04.anno/SNP.vcf"
    log:
        "logs/filter_maf_gene.log"
    shell:
        """
        module load plink/1.9
        plink --vcf {input} --maf 0.05 --geno 0.2 --recode vcf-iid --out result/04.anno/SNP > {log} 2>&1
        """


# 第一步：计算 LD 矩阵，生成保留的 SNP 列表
rule ld_indep_pairwise:
    input:
        "result/04.anno/SNP.vcf"
    output:
        "result/05.plink/ld.prune.in",
        "result/05.plink/ld.prune.out"
    log:
        "logs/ld_indep_pairwise.log"
    shell:
        """
        module load plink/1.9
        plink --vcf {input} --indep-pairwise 100 10 0.2 --out result/05.plink/ld > {log} 2>&1
        """

# 第二步：根据 LD 结果生成二进制文件格式
rule make_bed:
    input:
        vcf="result/04.anno/SNP.vcf",
        prune_in="result/05.plink/ld.prune.in"
    output:
        bed="result/05.plink/ld.prune.in.bed",
        bim="result/05.plink/ld.prune.in.bim",
        fam="result/05.plink/ld.prune.in.fam"
    log:
        "logs/make_bed.log"
    shell:
        """
        module load plink/1.9
        plink --vcf {input.vcf} --make-bed --extract {input.prune_in} --out result/05.plink/ld.prune.in > {log} 2>&1
        """

