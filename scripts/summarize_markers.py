import pandas as pd
from Bio import SeqIO
import os

# env: utils.yaml


# gather all the summary files and put them in a single table
summary_files = [file for file in snakemake.input if file.endswith(".summary")]
rows = list()
for summary_file in summary_files:
    lines = [line.strip().split("\t") for line in open(summary_file).readlines()]
    rows.append(lines[-1])

columns = ["contig", "TerL", "MCP", "portal", "primase"]

df = pd.DataFrame(rows, columns=columns)
df = df.sort_values(by="contig")

df.to_csv(snakemake.output.summary, sep="\t", index=False)


# store all the .faa files in the all_records dictionary
faa_files = [file for file in snakemake.input if file.endswith(".faa") and os.path.getsize(file) > 0]
all_records = dict()
for faa_file in faa_files:
    records = SeqIO.parse(faa_file, "fasta")
    for record in records:
        all_records[record.id] = record


# create the output dir for the .faa files
os.makedirs(snakemake.output.faa_dir, exist_ok = False)

# iterate the columns of the summary file while extracting the sequences
for column in columns[1:]: # don't look at the "contigs" column
    cells = df[column].tolist()
    # check if there were sequences detected for that marker. If so, add them to the final list
    to_write = [all_records[seq_id] for seq_id in cells if seq_id in all_records]
    if to_write:
        with open(f"{snakemake.output.faa_dir}/{column}.faa", "w") as fout:
            SeqIO.write(to_write, fout, "fasta")
