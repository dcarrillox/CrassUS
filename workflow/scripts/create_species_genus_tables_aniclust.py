from Bio import SeqIO
import pandas as pd

# env: compare_genomes.yaml


# -------------------------------------
# read aniclust results for genus level
lines = [line.strip().split("\t") for line in open(snakemake.input.gen[0])]
# to assign ids to the genera (genus__1, genus__2 ...), sort them by number of
# contigs in the cluster and assign lower ids to the clusters with more contigs
gen_cluster_counts = {line[0]:len(line[1].split(",")) for line in lines}
# sort the dict descending and assign a "genus__X" id to the repr genome
gen_cluster_counts = {k: v for k, v in sorted(gen_cluster_counts.items(), key=lambda item: item[1], reverse=True)}
cont = 1
repr_genus_ids = dict()
for repr in gen_cluster_counts:
    repr_genus_ids[repr] = f"genus__{cont}"
    cont += 1
# now assign the id to every genome in the genus cluster
genomes_genus_id = dict()
for line in lines:
    genomes = line[1].split(",")
    for genome in genomes:
        genomes_genus_id[genome] = repr_genus_ids[line[0]]



# -------------------------------------
# read aniclust results for species level.
lines = [line.strip().split("\t") for line in open(snakemake.input.sp[0])]

sp_cluster_counts = {line[0]:len(line[1].split(",")) for line in lines}
# sort the dict descending and assign a "sp__X" id to the repr genome
sp_cluster_counts = {k: v for k, v in sorted(sp_cluster_counts.items(), key=lambda item: item[1], reverse=True)}
cont = 1
repr_sp_ids = dict()
for repr in sp_cluster_counts:
    repr_sp_ids[repr] = f"species__{cont}"
    cont += 1
# now assign the id to every genome in the genus cluster.
# NB only repr genome of the species has genus id. Propagate the genus id to the
# rest of the species cluster
genomes_species_id = dict()
for line in lines:
    # get genus of the sp repr
    sp_repr_genus = genomes_genus_id[line[0]]

    genomes = line[1].split(",")
    for genome in genomes:
        genomes_species_id[genome] = repr_sp_ids[line[0]]
        genomes_genus_id[genome] = sp_repr_genus



# ---------------------------
# put in df and write to file
# read the fasta file with all the genomes, get their ids
all_genomes_ids = [record.id for record in SeqIO.parse(snakemake.input.fasta_all, "fasta")]

to_write = list()
for genome in all_genomes_ids:
    to_write.append([genome,
                    genomes_genus_id[genome],
                    genomes_species_id[genome]])

to_write_df = pd.DataFrame(to_write, columns=["genome", "genus", "species"])
to_write_df.set_index("genome", inplace=True)
# sort dataframe
to_write_df["ngen"] = to_write_df['genus'].str.split('__').str[-1].astype(int)
to_write_df["nsp"] = to_write_df['species'].str.split('__').str[-1].astype(int)
to_write_df.sort_values(by=["ngen", "nsp"], ascending=True, inplace=True)
to_write_df.drop(["ngen", "nsp"], axis = 1, inplace=True)

to_write_df.to_csv(snakemake.output[0], index=True, sep="\t")
