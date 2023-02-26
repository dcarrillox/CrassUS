from Bio import SearchIO, SeqIO
import pandas as pd

# start by reading the .faa file to get the the protein identifiers.
# store them along its n position in the genome
n_protein_ids = {record.id.split('|')[-1]:record.id for record in SeqIO.parse(snakemake.input.faa, "fasta")}

# read clusters_prots and assign a cluster_id to proteins, if they were part of a cluster
clusters_prots = {line.strip().split("\t")[1]:line.split("\t")[0] for line in open(snakemake.input.prots_clusters).readlines()[1:]}


# parse annotation with Yutin profiles
records = SearchIO.parse(snakemake.input.hmmtxt, "hmmer3-text")
profiles = {line.split("\t")[0].upper():line.strip().split("\t")[1:] for line in open(snakemake.params.yutin_names).readlines()}
#print(profiles)
yutin_annot = dict()

# init the dictionary for the annotation
for record in records:
    for hit in record.hits:
        yutin_annot[hit.id] = list()

# fill the dictionary in
records = SearchIO.parse(snakemake.input.hmmtxt, "hmmer3-text")
for record in records:
    for hit in record.hits:
        if hit.is_included and hit.bitscore > 15:
            #print([profiles[record.id.upper()], hit.evalue])
            yutin_annot[hit.id].append([profiles[record.id.upper()], hit.evalue])

# go through the hits and get the best full_length profile, or domain if there are no more hits
yutin_annot_final = dict()
for protein, hits in yutin_annot.items():
    if hits:
        sorted_hits = sorted(hits, key=lambda hit: float(hit[-1]))
        if len(sorted_hits) > 1:
            for hit in sorted_hits:
                if hit[0][1] == "full_length":
                    yutin_annot_final[protein] = hit[0][0]
                    break
        else:
            yutin_annot_final[protein] = sorted_hits[0][0][0]


# read gff file and create the genome table on the run
final_table = [["protein_id", "genome", "start", "end", "strand", "partial", "yutin", "cluster"]]

lines = [line.strip().split("\t") for line in open(snakemake.input.gff) if not line.startswith("#")]
for line in lines:
    n_gene = line[-1].split("ID=1_")[1].split(";")[0]
    # get start, stop, strand
    start,stop, strand = line[3], line[4], line[6]
    # check if there is annotation
    protein_id = n_protein_ids[n_gene]
    if protein_id in yutin_annot_final:
        yutin = yutin_annot_final[protein_id]
    else:
        yutin = ""
    # chekc if it is partial
    if ";partial=00;" in line[-1]:
        partial = "False"
    else:
        partial = "True"


    # check if the protein is in a cluster
    cluster_id = ""
    if protein_id in clusters_prots:
        cluster_id = clusters_prots[protein_id]

    final_table.append([protein_id, protein_id.split("|")[0], start, stop, strand, partial, yutin, cluster_id])


# write to output file
with open(snakemake.output[0], "w") as fout:
    for line in final_table:
        fout.write('\t'.join(line) + '\n')
