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

# don't remove query genome from the comparison list so there will be always an output file
#refs.remove(query)


if len(refs) > 3:
    refs = refs[:3]
# create execution (tmp) folder
os.makedirs(snakemake.params.tmp, exist_ok=True)
# copy fasta file to the execution folder
#os.system(f"cp {query} {' '.join(refs)} {snakemake.params.tmp}")
os.system(f"cp {' '.join(refs)} {snakemake.params.tmp}")

# run pyani
os.system(f"average_nucleotide_identity.py -i {snakemake.params.tmp} -o {snakemake.params.outdir} --workers {snakemake.threads} -m ANIb -f -v 2> {snakemake.log}") #
shutil.rmtree(snakemake.params.tmp, ignore_errors=True)
