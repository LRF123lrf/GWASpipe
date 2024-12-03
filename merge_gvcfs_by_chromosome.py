configfile: "config.yaml"

# 获取参考基因组和染色体列表
CHROMOSOMES = config["chromosomes"]

# 定义其他变量
REF = "ref/AGIS-1.0.fasta"
PFX = "AGIS-1.0"
RAWDATA = "data/samples"

(SAMPLES, READS) = glob_wildcards(RAWDATA + "/{sample}_{read}.fq.gz")
SAMPLES = sorted(set(SAMPLES))  # 去重并排序
print(SAMPLES)


rule all:
    input:
        expand("result/03.snp/{chrom}.merged.g.vcf.gz", chrom=CHROMOSOMES)

rule merge_gvcfs_by_chromosome:
    input:
        vcfs=expand("result/03.snp/{sample}_{chrom}.g.vcf",sample=SAMPLES, chrom=CHROMOSOMES),
        ref=REF
    output:
        "result/03.snp/{chrom}.merged.g.vcf.gz"
    params:
        vcf_list=lambda wildcards, input: " ".join(f"--variant {vcf}" for vcf in input.vcfs),
    shell:
        """
        module load GATK/4.6.0.0
        gatk CombineGVCFs \
            -R {input.ref} \
            {params.vcf_list} \
            -O {output}
        """

