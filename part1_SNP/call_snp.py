# Snakefile

# 加载配置文件
configfile: "config.yaml"

# 获取参考基因组和染色体列表
REF_GENOME = config["ref_genome"]
SAMPLES = config["samples"]
CHROMOSOMES = config["chromosomes"]

# 规则 all: 最终目标文件
rule all:
    input:
        expand("result/03.snp/{sample}_{chrom}.g.vcf", sample=SAMPLES, chrom=CHROMOSOMES)

# 规则1: 分染色体调用SNP (使用GATK)
rule call_variants_by_chromosome:
    input:
        bam="result/02.align/{sample}_sort_dup.bam",
        bai="result/02.align/{sample}_sort_dup.bam.bai",
        ref=REF_GENOME
    output:
        vcf="result/03.snp/{sample}_{chrom}.g.vcf"
    params:
        chrom="{chrom}"
    shell:
        """
        module load GATK/4.6.0.0
        gatk HaplotypeCaller -ERC GVCF \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            -L {params.chrom}
        """

