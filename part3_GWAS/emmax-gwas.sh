#BSUB -J pheno_kq.sh
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -o pheno_kq.sh.%J.out
#BSUB -e pheno_kq.sh.%J.err
#BSUB -q normal

/public/home/rfluo/01.software/emmax-intel64 -v -d 10 -t emmax_in -o output/pheno.kq -p Pheno.txt -k emmax_in.BN.kinf -c emmax.cov.txt
