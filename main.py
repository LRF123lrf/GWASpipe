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


# 定义所有规则
rule all:
    input:
        expand("{genome}.bwt", genome=REF),  # bwa索引
        expand("{genome}.fai", genome=REF),  # samtools索引
        expand("result/03.snp/{sample}_{chrom}.g.vcf", sample=SAMPLES, chrom=CHROMOSOMES),
        expand("result/03.snp/{chrom}.merged.g.vcf.gz", chrom=CHROMOSOMES),
        expand("result/05.plink/logs/log{K}.out", K=range(1, 9)),
        "ref/" + PFX + ".dict",  # GATK索引
        "result/03.snp/raw_variants.vcf.gz",  # 初步的合并为一个变异文件
        "result/03.snp/snp.vcf.gz",  # 只保留snp变异的文件
        "result/04.anno/anno.vcf.gz", #硬过滤、-q 0.1:major过滤、分染色体填充注释得到的文件
        "result/04.anno/SNP.vcf", #过滤--maf 0.05 --geno 0.2
        "input/SNP.vcf",
        "result/05.plink/ld.prune.in.bed",
        "result/05.plink/ld.prune.in.bim",
        "result/05.plink/ld.prune.in.fam",
        "result/06.pca/pca.bed",
        "result/06.pca/pca.bim",
        "result/06.pca/pca.fam",
        "result/06.pca/pca_output.eigenval",
        "result/06.pca/pca_output.eigenvec",
        "result/06.pca/pca_output.log",
        "result/06.pca/pca_output.nosex"


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
    params:
        vcf_list=lambda wildcards, input: " ".join(f"-I {vcf}" for vcf in input.genotyped_chr_vcfs)
    shell:
        """
        module load GATK/4.6.0.0
        gatk GatherVcfs \
            {params.vcf_list} \
            -O {output}
        """

rule index_vcf:
    input:
        vcf="result/03.snp/raw_variants.vcf.gz"
    output:
        idx="result/03.snp/raw_variants.vcf.gz.tbi"
    shell:
        """
        module load BCFtools/1.8
        bcftools index -t {input}
        """

rule gatk_SelectVariants:
    input:
        sample="result/03.snp/raw_variants.vcf.gz",
        index="result/03.snp/raw_variants.vcf.gz.tbi"
    output:
        snp="result/03.snp/snp.vcf.gz"
    params:
        genome=REF
    shell:
        """
        module load GATK/4.6.0.0
        gatk SelectVariants -R {params.genome} -V {input.sample} -select-type SNP -O {output.snp}
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

rule filter_pass_variants:
    input:
        "result/03.snp/snp.filtered.vcf.gz"
    output:
        "result/03.snp/snp.pass.vcf.gz"
    shell:
        """
        module load BCFtools/1.8
        bcftools view -i 'FILTER="PASS"' -Oz -o {output} {input}
        """

rule frequency_filter:
    input:
        "result/03.snp/snp.pass.vcf.gz"
    output:
        "result/04.anno/snp.maf.vcf.gz"
    threads: 10
    shell:
        """
        module load BCFtools/1.8
        bcftools norm -m -any {input} | \
        bcftools view -m 2 -M 2 -q 0.1:major -O z -o {output} --threads {threads}
        """

rule index_snp_maf_vcf:
    input:
        "result/04.anno/snp.maf.vcf.gz"
    output:
        "result/04.anno/snp.maf.vcf.gz.tbi"
    shell:
        """
        module load BCFtools/1.8
        bcftools index -t {input}
        """

rule split_by_chromosome:
    input:
        "result/04.anno/snp.maf.vcf.gz",
        tbi="result/04.anno/snp.maf.vcf.gz.tbi"
    output:
        "result/04.anno/split_by_chromosome/{chrom}.vcf.gz"
    params:
        chrom=lambda wildcards: wildcards.chrom
    log:
        "logs/split_by_chromosome/{chrom}.log"
    shell:
        """
        module load BCFtools/1.8
        bcftools view -r {params.chrom} {input} -O z -o {output} > {log} 2>&1
        """

rule genotype_imputation:
    input:
        "result/04.anno/split_by_chromosome/{chrom}.vcf.gz"
    output:
        "result/04.anno/{chrom}_imputed.vcf.gz"
    log:
        "logs/genotype_imputation/{chrom}.log"
    threads: 4
    params:
        out_prefix=lambda wildcards: f"result/04.anno/{wildcards.chrom}_imputed"
    shell:
        """
        module load beagle/4.1-Java-1.8.0_92
        java -jar $EBROOTBEAGLE/beagle.jar gt={input} impute=TRUE nthreads={threads} out={params.out_prefix} > {log} 2>&1
        """

rule index_imputed:
    input:
        "result/04.anno/{chrom}_imputed.vcf.gz"
    output:
        "result/04.anno/{chrom}_imputed.vcf.gz.tbi"
    log:
        "logs/index_imputed/{chrom}.log"
    shell:
        """
        module load BCFtools/1.8
        bcftools index -t {input} > {log} 2>&1
        """

rule annotate_variants:
    input:
        "result/04.anno/{chrom}_imputed.vcf.gz",
        tbi="result/04.anno/{chrom}_imputed.vcf.gz.tbi"
    output:
        "result/04.anno/{chrom}_anno.vcf.gz"
    log:
        "logs/annotate_variants/{chrom}.log"
    shell:
        """
        module load BCFtools/1.8
        bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {input} -Oz -o {output} > {log} 2>&1
        """

rule merge_annotated_chromosomes:
    input:
        expand("result/04.anno/{chrom}_anno.vcf.gz", chrom=CHROMOSOMES)
    output:
        "result/04.anno/anno.vcf.gz"
    log:
        "logs/merge_annotated_chromosomes.log"
    shell:
        """
        module load BCFtools/1.8
        bcftools concat -Oz -o {output} {input} > {log} 2>&1
        """

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


rule admixture:
    input:
        bed="result/05.plink/ld.prune.in.bed"
    output:
        "result/05.plink/logs/log{K}.out"
    shell:
        """
        /public/home/rfluo/01.software/admixture --cv {input.bed} {wildcards.K} | tee {output}
        mv ld* result/05.plink/logs/
        """


rule input_vcf:
    input:
        "result/04.anno/SNP.vcf"
    output:
        "input/SNP.vcf"
    log:
        "logs/input_vcf.log"
    shell:
        """
        ln -s ../{input} {output} > {log} 2>&1
        """


# 第一步：SNP.vcf文件转化格式
rule pca:
    input:
        "result/04.anno/SNP.vcf"
    output:
        bed="result/06.pca/pca.bed",
        bim="result/06.pca/pca.bim",
        fam="result/06.pca/pca.fam"
    log:
        "logs/pca.log"
    shell:
        """
        module load plink/1.9
        plink --vcf {input} --make-bed --out result/06.pca/pca > {log} 2>&1
        """


# 第二步：主成分分析
rule pca_analysis:
    input:
        bed="result/06.pca/pca.bed",
        bim="result/06.pca/pca.bim",
        fam="result/06.pca/pca.fam"
    output:
        eigenval="result/06.pca/pca_output.eigenval",
        eigenvec="result/06.pca/pca_output.eigenvec",
        log="result/06.pca/pca_output.log",
        nosex="result/06.pca/pca_output.nosex"
    log:
        "logs/pca_analysis.log"
    shell:
        """
        module load plink/1.9
        plink --bfile result/06.pca/pca --pca 30 --out result/06.pca/pca_output > {log} 2>&1
        """

