snakemake --jobs 999 --cluster-config cluster.yaml --cluster "sbatch -n {cluster.n}  -t {cluster.time}" --restart-times 2
