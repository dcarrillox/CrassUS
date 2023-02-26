#!/usr/bin/env python
# env: utils.yaml


import os
import pandas as pd
from functools import partial
import multiprocessing



def compute_genome_shared(all_genomes, clusters_list, n_prots_list, qgenome):
    # initialize the list for the genome that will contain the shared values
    shared_content = list()
    # iterate the genomes
    for tgenome in all_genomes:
        shared = 0
        if qgenome != tgenome:
            # iterate the clusters list of the query genome
            for i in range(0, len(clusters_list[qgenome])):
                # if both of the genomes have that ortholog...
                if clusters_list[qgenome][i] != 0 and clusters_list[tgenome][i] != 0:
                    # ... sum the n_prots of the query genome to shared
                    shared += int(clusters_list[qgenome][i])

            # # lowest n_protein of the two genomes in the denominator
            shared_perc = round(float(shared/n_prots_list[qgenome]),3)

            shared_content.append(shared_perc)
        # put 0 when self comparison
        else:
            shared_content.append(0)

    shared_content.insert(0, qgenome)
    return shared_content



# -------------------------------------------------------
# Assign an ID to each cluster with at least two proteins

# sort the tsv files by number of proteins in the cluster
os.system(f"cut -f1 {snakemake.input.tsv} | sort | uniq -c | sort -r -n > {snakemake.output.nprots}")

# create an identifier for each cluster, from larger to smaller cluster
cluster_ids = dict()
cont = 0
lines = [line.strip().split(" ") for line in open(snakemake.output.nprots).readlines()]
for line in lines:
    if int(line[0]) > 1:
        cont += 1
        cluster_ids[line[1]] = f"cl_{cont}"

# write new "table_clustering_ids.tsv" file with the cluster ids
cluster_ids_set = set(list(cluster_ids.keys()))
prots_clusters = [line.strip().split("\t") for line in open(snakemake.input.tsv).readlines() if line.split("\t")[0] in cluster_ids_set]
to_write = [[cluster_ids[prot[0]], prot[1]] for prot in prots_clusters]

to_write_df = pd.DataFrame(to_write, columns=["cluster_id", "protein"]).set_index("cluster_id")
to_write_df.to_csv(snakemake.output.table_clustering_ids, sep="\t")



# --------------------------------------------
# Create the presence-absence (presabs) matrix

# get all the genomes
prots = [line.strip().split("\t")[1] for line in open(snakemake.input.tsv).readlines()]
all_genomes = sorted(list(set([prot.split("|")[0] for prot in prots])))

# init an empty dataframe, rows are the protein clusters and columns the genomes
df = pd.DataFrame(0, index=[f"cl_{i}" for i in range(1, cont+1)], columns=all_genomes)

# iterate the tsv file and fill in the df adding 1 to the cluster-genome pair cell
lines = [line.strip().split("\t") for line in open(snakemake.input.tsv).readlines()]
# keep only clusters with more than one prot
lines = [line for line in lines if line[0] in cluster_ids]
for line in lines:
    df.loc[cluster_ids[line[0]], line[1].split("|")[0]] += 1

# write the presabs matrix to outfile
df.to_csv(snakemake.output.presabs, sep="\t")




# -------------------------------
# Calculate shared content matrix

# First get total n_prots per genome
genomes_totalp = {genome:0 for genome in all_genomes}
for prot in prots:
    genome_id = prot.split("|")[0]
    genomes_totalp[genome_id] += 1

# store each genome (column) in a dict, k=genome  v=clusters_presabs
genomes_clusters = {genome:df[genome].tolist() for genome in all_genomes}

# calculate sharing for all the genomes, both candidates and references
part = partial(compute_genome_shared, all_genomes, genomes_clusters, genomes_totalp)
pool = multiprocessing.Pool(processes=snakemake.threads)
shared_content = pool.map(part, all_genomes) # list of lists with the shared fraction
pool.close()
pool.join()

# construct df2 with the lists in shared_content
header = ["target_genome"] + all_genomes

df2 = pd.DataFrame(shared_content, columns=header).set_index("target_genome")
# transpose to have the contigs in the columns
df2 = df2.transpose()


# just before writing, add taxonomy of the reference
# read reference taxonomy
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()[1:]]
for line in lines:
    crass_taxonomy[line[0]] = {"family":line[1], "subfamily":line[2], "genus":line[3], "species":line[4]}

candidate_contigs = list()
for tgenome in df2.index:
    if tgenome in crass_taxonomy:
        df2.loc[tgenome, "family"] = crass_taxonomy[tgenome]["family"]
        df2.loc[tgenome, "subfamily"] = crass_taxonomy[tgenome]["subfamily"]
        df2.loc[tgenome, "genus"] = crass_taxonomy[tgenome]["genus"]
        df2.loc[tgenome, "species"] = crass_taxonomy[tgenome]["species"]
    else:
        candidate_contigs.append(tgenome)

# set the taxonomy columns as the first ones
taxa_columns = ["family", "subfamily", "genus", "species"]
df2 = df2[taxa_columns + [col for col in df2.columns if col not in taxa_columns]]

# write to final "shared_content_matrix_all.txt"
df2.to_csv(snakemake.output.shared_all, sep="\t")

# keep only candidate_contigs in the columns and write to "shared_content_matrix.txt"
df3 = df2[taxa_columns + sorted(candidate_contigs)]
df3.to_csv(snakemake.output.shared, sep="\t")
