#BSUB -J kinf.sh
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -o kinf.sh.%J.out
#BSUB -e kinf.sh.%J.err
#BSUB -q normal

module load plink/1.9

#将基因型文件转换为EMMAX格式，注意染色体数目的设置。
plink --vcf SNP.vcf --recode 12 transpose --output-missing-genotype 0 --out emmax_in --chr-set 18 --allow-extra-chr
#得到tfam文件

#群体亲缘关系分析
/public/home/rfluo/01.software/emmax-kin-intel64 emmax_in -v -d 10 -o emmax_in.BN.kinf
