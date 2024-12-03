rule all:
    input:
        "result/05.plink/logs/admixture_plot.png"  # 最终目标文件

rule plot_admixture:
    input:
        Q_file="result/05.plink/logs/ld.prune.in.3.Q"  # 输入的 Q 文件
    output:
        "result/05.plink/logs/admixture_plot.png"  # 输出的图像文件
    log:
        "logs/plot_admixture.log"
    shell:
        """
        module load Python/3.6.5
        python plot_admixture.py {input.Q_file}
        """

