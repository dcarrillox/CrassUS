from Bio import SeqIO
import os

# env: utils.yaml

os.makedirs(snakemake.output[0], exist_ok = True)

with open(snakemake.input[0]) as fin:
    contigs = [line.strip() for line in fin.readlines()]

# get the different samples for which contigs were detected
samples = list(set([contig.split("_")[0] for contig in contigs]))

# iterate the samples
for sample in samples:
    # read the assembly file for the sample
    assembly = f"{snakemake.params.assemblies_dir}/{sample}.fasta"
    records = {record.id:record for record in SeqIO.parse(assembly, "fasta")}

    # now iterate the contigs. For those coming from the sample, grab their sequence
    for contig in contigs:
        sample_contig = contig.split("_")[0]

        if sample_contig == sample:
            out_file = f"{snakemake.params.contigs_dir}/{contig}.fasta"
            with open(out_file, "w") as fout:
                SeqIO.write(records[contig], fout, "fasta")


## update genomes list for the fastANI step later.
# first read the list with the genomes from the database
lines = [line.strip() for line in open(snakemake.params.genomes_list).readlines()]
# then update the list with the found genomes
for contig in contigs:
    lines.append(f"{snakemake.params.contigs_dir}/{contig}.fasta")
# write to file
with open(snakemake.output.fastani_list, "w") as fout:
    for line in lines:
        fout.write(f"{line}\n")
