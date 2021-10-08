import os
import pandas as pd
import shutil


# read with pandas while sorting by cov
df = pd.read_csv(snakemake.input.ani[0], sep="\t", header=None,
                names=["query", "ref", "ani", "cov", "total"]).sort_values(
                ["cov", "ani"], ascending=[False, False], ignore_index=True)

# get query and top 5 covered refs
query = df["query"][0]
refs   = df["ref"].to_list()

if len(refs) > 5:
    refs = refs[:5]

# check that query contig is one of the top refs. This is not always the case, specially
# at short contigs that are highly covered by itself but also by other longer contigs
if query not in refs:
    refs.insert(0,query)


# create execution (tmp) folder
os.makedirs(snakemake.params.tmp, exist_ok=True)
# copy fasta file to the execution folder
#os.system(f"cp {query} {' '.join(refs)} {snakemake.params.tmp}")
os.system(f"cp {' '.join(refs)} {snakemake.params.tmp}")

# run pyani
os.system(f"average_nucleotide_identity.py -i {snakemake.params.tmp} -o {snakemake.params.outdir} --workers {snakemake.threads} -m ANIb -f -v 2> {snakemake.log}") #
shutil.rmtree(snakemake.params.tmp, ignore_errors=True)
