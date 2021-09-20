import os
import pandas as pd
import shutil


# read with pandas while sorting by cov
df = pd.read_csv(snakemake.input.ani[0], sep="\t", header=None,
                names=["query", "ref", "ani", "cov", "total"]).sort_values(
                "cov", ascending=False, ignore_index=True)
# get query and top 5 covered refs
query = df["query"][0]
refs   = df["ref"].to_list()

refs.remove(query)


# if there was any hit
if refs:
    if len(refs) > 3:
        refs = refs[:3]
    # create execution (tmp) folder
    os.makedirs(snakemake.params.tmp, exist_ok=True)
    # copy fasta file to the execution folder
    os.system(f"cp {query} {' '.join(refs)} {snakemake.params.tmp}")

    # run pyani
    os.system(f"average_nucleotide_identity.py -i {snakemake.params.tmp} -o {snakemake.params.outdir} --workers {snakemake.threads} -m ANIb -f -v 2> {snakemake.log}") #
    os.system(f"touch {snakemake.output.done}")
    shutil.rmtree(snakemake.params.tmp, ignore_errors=True)

else:
    os.makedirs(snakemake.params.outdir, exist_ok=True)
    os.system(f"touch {snakemake.output.done}")
