import os
import pandas as pd


# get the query contig
query_contig = snakemake.wildcards.contig

# create output directory

os.makedirs(snakemake.output.tmp, exist_ok=True)

# check if there are pyani coverage results
coverage_file = snakemake.input[0]

# read coverage file
df = pd.read_csv(coverage_file, sep="\t", index_col=0)
tdf = df.transpose()
# sort df by the query_contig column to have them sorted
tdf = tdf.sort_values(query_contig,  ascending=False)
# get sorted coverage values for the query contig
coverages_sorted = tdf[query_contig]
# most covering genome should be the genome itself, second should be the chosen
# for the comparison. But NW some genomes might not be similar to anything else
# but to themselves, so length of the list might be only 1.
if len(coverages_sorted) > 1:
    target_genome = coverages_sorted.index.tolist()[1]
else:
    target_genome = coverages_sorted.index.tolist()[0]

# genome can be from the reference or from the found genomes. Check if the files
# for both scenarios are present
if os.path.isfile(f"{snakemake.params.refgenomes_dir}/{target_genome}.fasta"):
    target_genome_fasta = f"{snakemake.params.refgenomes_dir}/{target_genome}.fasta"
else:
    target_genome_fasta = f"{snakemake.params.contigs_dir}/{target_genome}.fasta"
# copy target genome to tblastx results folder and make a blast db of it
os.system(f"cp {target_genome_fasta} {snakemake.output.tmp}")
os.system(f"makeblastdb -in {snakemake.output.tmp}/{target_genome}.fasta -out {snakemake.output.tmp}/{target_genome} -dbtype nucl")

# run megablast
query_contig_fasta = f"{snakemake.params.contigs_dir}/{query_contig}.fasta"
db = f"{snakemake.output.tmp}/{target_genome}"
outfile = f"{snakemake.output.tmp}/{query_contig}.megablast"
os.system(f"blastn -query {query_contig_fasta} -db {db} -out {outfile} "
"-task megablast -dust no -soft_masking false "
"-outfmt '6 qseqid sseqid pident length mismatch nident qlen qstart qend slen sstart send evalue bitscore'")
# move file to its final location
os.system(f"mv {outfile} {snakemake.output.megablast}")
