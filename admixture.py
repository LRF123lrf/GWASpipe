rule all:
    input:
        expand("result/05.plink/logs/log{K}.out", K=range(1, 9))

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

