from Bio import SeqIO, SearchIO

# Looking at hit's evalue might be too greedy, better check the ivalue of the hsp(s)
records = SearchIO.parse(snakemake.input.hmmtxt, "hmmer3-text")

terL_profiles = ["VP00412", "VP02579", "VP02686"]

# identify hits to the terminase profiles
terL_hits = list()
for record in records:
    if record.id in terL_profiles:
        for hit in record.hits:
            for hsp in hit.hsps:
                if hsp.is_included:
                    terL_hits.append(hit.id)

# check if any terL hits were identified
to_write_summary = list()

if terL_hits:
    # make the hits uniq
    terL_hits = list(set(terL_hits))
    # init a list to store the ids for the faa sequence
    terL_faa  = list()

    # read gff file
    gff_lines = [line.strip() for line in open(snakemake.input.gff).readlines() if not line.startswith("#")]

    # check that line gene_id - 1 in the gff is not partial
    to_write_summary_summary = list()
    for gene in terL_hits:
        n_gene = int(gene.split("|")[-1])
        if "partial=00" in gff_lines[n_gene-1]:
            to_write_summary.append(f"{gene}\tcomplete TerL")
            terL_faa.append(gene)
        else:
            to_write_summary.append(f"{gene}\tpartial TerL")

else:
    to_write_summary("No TerL hits")


# write summary file
with open(snakemake.output.summary, "w") as fout:
    for line in to_write_summary:
        fout.write(f"{line}\n")


# write .faa file
with open(snakemake.output.faa, "w") as fout:
    # if there is a TerL sequence(s) to grab from the prodigal file...
    if terL_faa:
        to_write_faa = list()
        records = {record.id:record for record in SeqIO.parse(snakemake.input.faa, "fasta")}
        for seq_id in terL_faa:
            to_write_faa.append(records[seq_id])
        SeqIO.write(to_write_faa, fout, "fasta")
