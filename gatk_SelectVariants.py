# 加载配置文件
configfile: "config.yaml"

# 获取参考基因组和染色体列表
CHROMOSOMES = config["chromosomes"]

# 定义其他变量
REF = "ref/AGIS-1.0.fasta"
PFX = "AGIS-1.0"
RAWDATA = "data"

(SAMPLES, READS) = glob_wildcards(RAWDATA + "/{sample}_{read}.fq.gz")
SAMPLES = sorted(set(SAMPLES))  # 去重并排序
print(SAMPLES)

rule all:
    input:
        "result/03.snp/snp.vcf.gz"

rule gatk_SelectVariants:
    input:
        sample="result/03.snp/raw_variants.vcf.gz"
    output:
        snp="result/03.snp/snp.vcf.gz"
    params:
        genome=REF
    shell:
        """
        module load GATK/4.6.0.0
        gatk SelectVariants -R {params.genome} -V {input.sample} -select-type SNP -O {output.snp}
        """
