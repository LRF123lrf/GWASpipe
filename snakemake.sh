module load snakemake/5.4.0
snakemake --rerun-incomplete --cluster-config cluster.json --cluster  \
'bsub -n {cluster.nCPUs} -q {cluster.queue} -R {cluster.resources}'  -j 200  -s main.py
