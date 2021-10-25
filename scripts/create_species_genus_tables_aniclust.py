from Bio import SeqIO
import pandas as pd

# env: compare_genomes.yaml

# read the fasta file with all the genomes, get their ids
all_genomes_ids = [record.id for record in SeqIO.parse(snakemake.input.fasta_all, "fasta")]

# read genus
gen_df = pd.read_csv(snakemake.input.gen[0], sep="\t", names=["repr", "cluster"], header=None,  index_col=0)

cont = 1
genomes_genus_id = dict()
for genome in gen_df.index:
    cluster = gen_df.loc[genome, "cluster"].split(",")
    for cluster_genome in cluster:
        genomes_genus_id[cluster_genome] = f"genus__{cont}"
    cont += 1


# read species
sp_df = pd.read_csv(snakemake.input.sp[0], sep="\t", names=["repr", "cluster"], header=None,  index_col=0)

cont = 1
genomes_species_id = dict()
genomes_order = list()
for genome in sp_df.index:
    cluster = sp_df.loc[genome, "cluster"].split(",")
    for cluster_genome in cluster:
        genomes_species_id[cluster_genome] = f"species__{cont}"
        genomes_genus_id[cluster_genome] = genomes_genus_id[genome]
        genomes_order.append(cluster_genome)
    cont += 1


# put in df and write to file
to_write = list()
for genome in all_genomes_ids:
    to_write.append([genome,
                    genomes_genus_id[genome],
                    genomes_species_id[genome]])


#
print(genomes_order)
to_write_df = pd.DataFrame(to_write, columns=["genome", "genus", "species"])
to_write_df.set_index("genome", inplace=True)
to_write_df = to_write_df.reindex(genomes_order)
to_write_df.to_csv(snakemake.output[0], index=True, sep="\t")
