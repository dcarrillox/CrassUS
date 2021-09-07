from Bio import SeqIO
import os

# env: utils.yaml

os.makedirs(snakemake.output[0], exist_ok = True)

with open(snakemake.input[0]) as fin:
    contigs = [line.strip() for line in fin.readlines()]


for contig in contigs:
    sample = contig.split("_")[0]

    # read the assembly file for the sample
    assembly = f"{snakemake.params.assemblies_dir}/{sample}.fasta"
    records = {record.id:record for record in SeqIO.parse(assembly, "fasta")}

    out_file = f"{snakemake.params.contigs_dir}/{contig}.fasta"
    with open(out_file, "w") as fout:
        SeqIO.write(records[contig], fout, "fasta")
