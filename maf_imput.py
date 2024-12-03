configfile: "config.yaml"
CHROMOSOMES = config["chromosomes"]


rule all:
    input:
        "result/04.anno/anno.vcf.gz"

rule filter_pass_variants:
    input:
        "result/03.snp/snp.filtered.vcf.gz"
    output:
        "result/03.snp/snp.pass.vcf.gz"
    log:
        "logs/filter_pass_variants.log"
    shell:
        """
        module load BCFtools/1.8
        bcftools view -i 'FILTER="PASS"' -Oz -o {output} {input} > {log} 2>&1
        """

rule frequency_filter:
    input:
        "result/03.snp/snp.pass.vcf.gz"
    output:
        "result/04.anno/snp.maf.vcf.gz"
    log:
        "logs/frequency_filter.log"
    threads: 10
    shell:
        """
        module load BCFtools/1.8
        bcftools norm -m -any {input} | \
        bcftools view -m 2 -M 2 -q 0.1:major -O z -o {output} --threads {threads} > {log} 2>&1
        """

rule index_vcf:
    input:
        "result/04.anno/snp.maf.vcf.gz"
    output:
        "result/04.anno/snp.maf.vcf.gz.tbi"
    log:
        "logs/index_vcf.log"
    shell:
        """
        module load BCFtools/1.8
        bcftools index -t {input} > {log} 2>&1
        """

rule split_by_chromosome:
    input:
        "result/04.anno/snp.maf.vcf.gz"
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
    threads: 6
    params:
        out_prefix=lambda wildcards: f"result/04.anno/{wildcards.chrom}_imputed"
    shell:
        """
        module load beagle/4.1-Java-1.8.0_92
        java -jar $EBROOTBEAGLE/beagle.jar gt={input} impute=TRUE nthreads={threads} out={params.out_prefix} > {log} 2>&1
        """

rule index_imputed:
    input:
        "result/04.anno/{chrom}.imputed.vcf.gz"
    output:
        "result/04.anno/{chrom}.imputed.vcf.gz.tbi"
    log:
        "logs/index_imputed/{chrom}.log"
    shell:
        """
        module load BCFtools/1.8
        bcftools index -t {input} > {log} 2>&1
        """

rule annotate_variants:
    input:
        "result/04.anno/{chrom}_imputed.vcf.gz"
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

