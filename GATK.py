# 加载配置文件
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


# 定义所有规则
rule all:
    input:
        expand("{genome}.bwt", genome=REF),  # bwa索引
        expand("{genome}.fai", genome=REF),  # samtools索引
        expand("result/03.snp/{sample}_{chrom}.g.vcf", sample=SAMPLES, chrom=CHROMOSOMES),
        expand("result/03.snp/{chrom}.merged.g.vcf.gz", chrom=CHROMOSOMES),
        "ref/" + PFX + ".dict",  # GATK索引
        "result/03.snp/raw_variants.vcf.gz",  # 最终的全基因组合并文件
        "result/03.snp/snp.vcf.gz"  # 最终的全基因组合并文件

rule fastp:
    input:
        R1=RAWDATA + "/{sample}_1.fq.gz",
        R2=RAWDATA + "/{sample}_2.fq.gz"
    output:
        R1="result/01.data/trim/{sample}/{sample}_1_val_1.fq.gz",
        R2="result/01.data/trim/{sample}/{sample}_2_val_2.fq.gz"
    params:
        dirs="result/01.data/trim/{sample}"
    log:
        "logs/trim/{sample}.log"
    shell:
        """
        module load TrimGalore/0.6.6
        trim_galore -j {threads} -q 25 --phred33 --length 50 -e 0.1 --stringency 3 --paired -o {params.dirs} --fastqc {input.R1} {input.R2} 
        """

rule bwa_index:
    input:
        genome=REF
    output:
        "{genome}.bwt"
    shell:
        """
        module load BWA/0.7.17
        bwa index {input.genome}
        """

rule bwa_map:
    input:
        R1="result/01.data/trim/{sample}/{sample}_1_val_1.fq.gz",
        R2="result/01.data/trim/{sample}/{sample}_2_val_2.fq.gz",
        genome_index=expand("{genome}.bwt", genome=REF)
    output:
        "result/02.align/{sample}_sort.bam"
    params:
        genome=REF
    shell:
        """
        module load BWA/0.7.17
        module load SAMtools/1.9
        bwa mem -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:Abace' -t {threads} {params.genome} {input.R1} {input.R2} | samtools sort -@ {threads} -o {output}
        """

rule samtools_faidx:
    input:
        genome=REF
    output:
        "{genome}.fai"
    shell:
        """
        module load SAMtools/1.9
        samtools faidx {input.genome}
        """

rule GATK_dict:
    input:
        genome=REF
    output:
        dict="data/ref/{genome_profix}.dict"
    shell:
        """
        module load GATK/4.6.0.0
        gatk CreateSequenceDictionary -R {input.genome} -O {output.dict}
        """

rule GATK_MarkDuplicates:
    input:
        sample="result/02.align/{sample}_sort.bam"
    output:
        dup_bam="result/02.align/{sample}_sort_dup.bam"
    log:
        metrics="result/02.align/log/{sample}_markdup_metrics.txt"
    shell:
        """
        module load GATK/4.6.0.0
        gatk MarkDuplicates -I {input.sample} -O {output.dup_bam} -M {log.metrics}
        """

rule samtools_index:
    input:
        "result/02.align/{sample}_sort_dup.bam"
    output:
        "result/02.align/{sample}_sort_dup.bam.bai"
    shell:
        """
        module load SAMtools/1.9
        samtools index {input}
        """

# 规则1: 分染色体调用SNP (使用GATK)
rule call_variants_by_chromosome:
    input:
        bam="result/02.align/{sample}_sort_dup.bam",
        bai="result/02.align/{sample}_sort_dup.bam.bai",
        ref=REF
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


# 规则2: 按染色体合并所有样本的.g.vcf文件，并压缩为.g.vcf.gz
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

# 规则3:转化各染色体.g.vcf.gz文件为.vcf.gz文件
rule GATK_GenotypeGVCFs:
    input:
        combine_chr_vcf="result/03.snp/{chrom}.merged.g.vcf.gz",
        ref=REF
    output:
        genotyped_chr_vcf="result/03.snp/{chrom}.merged.vcf.gz"
    params:
        chrom="{chrom}"
    shell:
        """
        module load GATK/4.6.0.0
        gatk GenotypeGVCFs -R {input.ref} -V {input.combine_chr_vcf} -O {output.genotyped_chr_vcf}
        """

# 规则4:合并为一个.vcf.gz文件
rule GATK_GatherVcfs:
    input:
        genotyped_chr_vcfs=expand("result/03.snp/{chrom}.merged.vcf.gz", chrom=CHROMOSOMES)
    output:
        "result/03.snp/raw_variants.vcf.gz"
    log:
        "logs/gather_vcfs.log"
    params:
        vcf_list=lambda wildcards, input: " ".join(f"-I {vcf}" for vcf in input.genotyped_chr_vcfs)
    shell:
        """
        module load GATK/4.6.0.0
        gatk GatherVcfs \
            {params.vcf_list} \
            -O {output} > {log} 2>&1
        """

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
        gatk SelectVariants -R {params.genome} -V {input.sample} -O {output.snp} --select-type-to-include SNP
        """

rule gatk_VariantFiltration:
    input:
        sample="result/03.snp/snp.vcf.gz"
    output:
        snp="result/03.snp/snp.filtered.vcf.gz"
    params:
        genome=REF
    shell:
        """
        module load GATK/4.6.0.0
        gatk VariantFiltration -R {params.genome} -V {input.sample} -O {output.snp} --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum< -12.5 || ReadPosRankSum < -8.0" --filter-name "filter"
        """

