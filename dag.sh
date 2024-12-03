snakemake --dag -s main.py > dag.dot
dot -Tpng dag.dot -o dag.png
