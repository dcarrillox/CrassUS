import os
import pandas as pd


# get the query contig
query_contig = snakemake.wildcards.contig

# create output directory

os.makedirs(snakemake.params.outdir, exist_ok=True)

# check if there are pyani coverage results
coverage_file = f"{snakemake.params.pyani_dir}/ANIb_alignment_coverage.tab"
if os.path.isfile(coverage_file):
    # read coverage file
    df = pd.read_csv(coverage_file, sep="\t", index_col=0)
    tdf = df.transpose()
    # sort df by the query_contig column to have them sorted
    tdf = tdf.sort_values(query_contig,  ascending=False)
    # get sorted coverage values for the query contig
    coverages_sorted = tdf[query_contig]
    # most covering genome should be the genome itself, second should be the chosen
    # for the comparison
    target_genome = coverages_sorted.index.tolist()[1]

    # genome can be from the reference or from the found genomes. Check if the files
    # for both scenarios are present
    if os.path.isfile(f"{snakemake.params.refgenomes_dir}/{target_genome}.fasta"):
        target_genome_fasta = f"{snakemake.params.refgenomes_dir}/{target_genome}.fasta"
    else:
        target_genome_fasta = f"{snakemake.params.contigs_dir}/{target_genome}.fasta"
    # copy target genome to tblastx results folder and make a blast db of it
    os.system(f"cp {target_genome_fasta} {snakemake.params.outdir}")
    os.system(f"makeblastdb -in {snakemake.params.outdir}/{target_genome}.fasta -out {snakemake.params.outdir}/{target_genome} -dbtype nucl")

    # run tblastx
    query_contig_fasta = f"{snakemake.params.contigs_dir}/{query_contig}.fasta"
    db = f"{snakemake.params.outdir}/{target_genome}"
    outfile = f"{snakemake.params.outdir}/{query_contig}.megablast"
    os.system(f"blastn -query {query_contig_fasta} -db {db} -out {outfile} "
    "-task megablast -dust no -soft_masking false "
    "-outfmt '6 qseqid sseqid pident length mismatch nident qlen qstart qend slen sstart send evalue bitscore'")

    # write the mock output file
    os.system(f"touch {snakemake.output.done}")


else:
    os.system(f"touch {snakemake.output.done}")
