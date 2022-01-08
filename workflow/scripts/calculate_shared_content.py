import os
import pandas as pd
from functools import partial
import multiprocessing



################
# presabs matrix
# presabs matrix is generated for ALL the genomes: ref + found
################

# get all the genomes
prots = [line.strip().split("\t")[1] for line in open(snakemake.input.tsv).readlines()]
all_genomes = sorted(list(set([prot.split("|")[0] for prot in prots])))

# sort the tsv files by number of proteins in the cluster
os.system(f"cut -f1 {snakemake.input.tsv} | sort | uniq -c | sort -r -n > {snakemake.output.nprots}")

# create an identifier for each cluster, from biggest cluster to smallest
cluster_ids = dict()
cont = 0
lines = [line.strip().split(" ") for line in open(snakemake.output.nprots).readlines()]
for line in lines:
    if int(line[0]) > 1:
        cont += 1
        cluster_ids[line[1]] = f"cluster_{cont}"

# write new table_clustering_ids.tsv file, which is like the original table_clustering.tsv
# file but with the representative replaced by the cluster_id. Another difference is that
# it only contains clusters with more than two sequences
print("Annotating 'table_clustering.tsv' with cluster_ids...")

cluster_ids_set = set(list(cluster_ids.keys()))
prots_clusters = [line.strip().split("\t") for line in open(snakemake.input.tsv).readlines() if line.split("\t")[0] in cluster_ids_set]
to_write = [[cluster_ids[prot[0]], prot[1]] for prot in prots_clusters]

to_write_df = pd.DataFrame(to_write, columns=["cluster_id", "protein"]).set_index("cluster_id")
to_write_df.to_csv(snakemake.output.table_clustering_ids, sep="\t")

print("done")


# init an empty dataframe, rows are the clusters and columns the genomes
df = pd.DataFrame(index=[f"cluster_{i}" for i in range(1, cont+1)], columns=all_genomes)
df = df.fillna(0)


# iterate the tsv file and fill in the df adding 1 to the cluster-genome pair cell
lines = [line.strip().split("\t") for line in open(snakemake.input.tsv).readlines()]
# keep only clusters with more than one prot
lines = [line for line in lines if line[0] in cluster_ids]
for line in lines:
    df.loc[cluster_ids[line[0]], line[1].split("|")[0]] += 1

# write the presabs matrix to outfile
df.to_csv(snakemake.output.presabs, sep="\t")


################
# shared content
################

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



# get total n_prots per genome
genomes_totalp = {genome:0 for genome in all_genomes}
for prot in prots:
    genome_id = prot.split("|")[0]
    genomes_totalp[genome_id] += 1

# store each genome (column) in a dict, k=genome  v=clusters_presabs
genomes_clusters = {genome:df[genome].tolist() for genome in all_genomes}

## TESTING MULTITHREADING
print(f"Calculating shared content with {snakemake.threads} CPUs...")
# prepare the partial function
part = partial(compute_genome_shared, all_genomes, genomes_clusters, genomes_totalp)

pool = multiprocessing.Pool(processes=snakemake.threads)
# calculate sharing for all the genomes instead, found + ref
#contigs = sorted([os.path.basename(prot).split("_tbl-")[0] for prot in snakemake.input.prots_files])
#shared_content = pool.map(part, contigs)
shared_content = pool.map(part, all_genomes)
pool.close()
pool.join()

# construct df2 with the lists in shared_content
header = ["target_genome"] + all_genomes

df2 = pd.DataFrame(shared_content, columns=header).set_index("target_genome")
# transpose to have the contigs in the columns
df2 = df2.transpose()
print(df2.shape)
print("done")

## BEFORE TESTING MULTITHREADING
# # init a dataframe, all genomes in the rows, found contigs in the columns
# contigs = sorted([os.path.basename(prot).split("_tbl-")[0] for prot in snakemake.input.prots_files])
# df2 = pd.DataFrame(index=all_genomes, columns=contigs)
# df2 = df2.fillna(0)
#
#
# # iterate the genomes and compare them
# for contig in contigs:
#     for ref in all_genomes:
#         shared = 0
#         # check it is not a self comparison
#         if ref != contig:
#             for i in range(0, cont):
#                 if genomes_clusters[contig][i] > 0 and genomes_clusters[ref][i] > 0:
#                     #print(genomes_clusters[contig][i], genomes_clusters[ref][i])
#                     shared += genomes_clusters[contig][i]
#
#             # compute the shared perc
#             shared_perc = round(float(shared/genomes_totalp[contig]), 3)
#             # fill df2 with it. query is the row, ref is the column
#             df2.loc[ref, contig] = shared_perc


# write the shared content matrix to outfile
# NW that only found contigs are in the columns. The table is not intended to be
# simetric, it only contains sharing values of the contigs but not of the reference
# genomes
# just before writing, add taxonomy of the reference
# read reference taxonomy
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()[1:]]
for line in lines:
    crass_taxonomy[line[0]] = {"family":line[1], "subfamily":line[2], "genus":line[3], "species":line[4]}

for tgenome in df2.index:
    if tgenome in crass_taxonomy:
        df2.loc[tgenome, "family"] = crass_taxonomy[tgenome]["family"]
        df2.loc[tgenome, "subfamily"] = crass_taxonomy[tgenome]["subfamily"]
        df2.loc[tgenome, "genus"] = crass_taxonomy[tgenome]["genus"]
        df2.loc[tgenome, "species"] = crass_taxonomy[tgenome]["species"]


taxa_columns = ["family", "subfamily", "genus", "species"]
df2 = df2[taxa_columns + [ col for col in df2.columns if col not in taxa_columns]]



df2.to_csv(snakemake.output.shared, sep="\t")
