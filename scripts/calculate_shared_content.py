import os
import pandas as pd

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

# get total n_prots per genome
genomes_totalp = {genome:0 for genome in all_genomes}
for prot in prots:
    genome_id = prot.split("|")[0]
    genomes_totalp[genome_id] += 1

# store each genome (column) in a dict, k=genome  v=clusters_presabs
genomes_clusters = {genome:df[genome].tolist() for genome in all_genomes}

# init a dataframe, all genomes in the rows, found contigs in the columns
contigs = sorted([os.path.basename(prot).split("_tbl-")[0] for prot in snakemake.input.prots_files])
df2 = pd.DataFrame(index=all_genomes, columns=contigs)
df2 = df2.fillna(0)


# iterate the genomes and compare them
for contig in contigs:
    for ref in all_genomes:
        shared = 0
        # check it is not a self comparison
        if ref != contig:
            for i in range(0, cont):
                if genomes_clusters[contig][i] > 0 and genomes_clusters[ref][i] > 0:
                    #print(genomes_clusters[contig][i], genomes_clusters[ref][i])
                    shared += genomes_clusters[contig][i]

            # compute the shared perc
            shared_perc = round(float(shared/genomes_totalp[contig]), 3)
            # fill df2 with it. query is the row, ref is the column
            df2.loc[ref, contig] = shared_perc


# write the shared content matrix to outfile
# NW that only found contigs are in the columns. The table is not intended to be
# simetric, it only contains sharing values of the contigs but not of the reference
# genomes
# just before writing, add taxonomy of the reference
# read reference taxonomy
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
for line in lines:
    crass_taxonomy[line[0]] = {"family":line[2], "subfamily":line[3], "genus":line[4]}

for tgenome in df2.index:
    if tgenome in crass_taxonomy:
        df2.loc[tgenome, "family"] = crass_taxonomy[tgenome]["family"]
        df2.loc[tgenome, "subfamily"] = crass_taxonomy[tgenome]["subfamily"]
        df2.loc[tgenome, "genus"] = crass_taxonomy[tgenome]["genus"]

taxa_columns = ["family", "subfamily", "genus"]
df2 = df2[taxa_columns + [ col for col in df2.columns if col not in taxa_columns]]



df2.to_csv(snakemake.output.shared, sep="\t")
