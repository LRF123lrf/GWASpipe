configfile: "config.yaml"
CHROMOSOMES = config["chromosomes"]
REF = "ref/AGIS-1.0.fasta"
PFX = "AGIS-1.0"
RAWDATA = "data"

(SAMPLES, READS) = glob_wildcards(RAWDATA + "/{sample}_{read}.fq.gz")
SAMPLES = sorted(set(SAMPLES))

rule all:
    input:
        expand("result/03.snp/{chrom}.merged.g.vcf.gz", chrom=CHROMOSOMES)

# 合并每条染色体的 g.vcf 文件
rule merge_gvcfs_by_chromosome:
    input:
        vcfs=lambda wildcards: expand("result/03.snp/{sample}_{chrom}.g.vcf", 
                                      sample=SAMPLES, chrom=wildcards.chrom),
        ref=REF
    output:
        "result/03.snp/{chrom}.merged.g.vcf.gz"
    params:
        vcf_list=lambda wildcards, input: " ".join(f"--variant {vcf}" for vcf in input.vcfs)
    log:
        "logs/merge_gvcfs/{chrom}.log"
    shell:
        """
        module load GATK/4.6.0.0
        gatk CombineGVCFs \
            -R {input.ref} \
            {params.vcf_list} \
            -O {output} > {log} 2>&1
        """
