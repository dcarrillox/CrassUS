import pandas as pd

# envs: compare_genomes.yaml

# ------------------------------------------------------------------------------
# read Crassvirales reference taxonomy, store genomes and their genus and species
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()[1:]]
crass_taxonomy = dict()
for line in lines:
    crass_taxonomy[line[0]] = {"genus":line[3], "species":line[4]}



# ---------------------
# read assignments file and store which genus__XX and species__XX represent
# reference taxa
lines = [line.strip().split("\t") for line in open(snakemake.input.assignments[0]).readlines()[1:]]
# Genus rank
genid_ref = {line[1]:list() for line in lines}
for line in lines:
    if line[0] in crass_taxonomy:
        genid_ref[line[1]] += [crass_taxonomy[line[0]]["genus"]]

# Species rank
spid_ref = {line[2]:list() for line in lines}
for line in lines:
    if line[0] in crass_taxonomy:
        spid_ref[line[2]] += [crass_taxonomy[line[0]]["species"]]




# read anicalc table so I can get later the most similar ref genome
anicalc_df = pd.read_csv(snakemake.input.anicalc[0], sep="\t", header=0, index_col=0)
anicalc_genomes = set(anicalc_df.index.tolist())



# ----
# iterate assignments and anicalc tables
# init a dict to store the final results. Include only candidates (so exclude
# reference Crassvirales)
final_dict = {line[0]:
                        {"genus_cluster":str(),
                        "genus_names":str(),
                        "species_cluster":str(),
                        "species_names":str(),
                        "most_similar_genome":str(),
                        "most_similar_ref_genome":str(),
                        "ref_genus":str(),
                        "ref_species":str(),
                        "pid":str(),
                        "qcov":0}
                for line in lines if line[0] not in crass_taxonomy and line[0] in anicalc_genomes}

# iterate assignments
for line in lines:
    qgenome = line[0]

    if line[0] in final_dict:
        # Genus
        final_dict[qgenome]["genus_cluster"] = line[1]
        genus_names = list(set(genid_ref[line[1]]))
        if genus_names:
            if len(genus_names) == 1:
                genus_name = genus_names[0]
            else:
                genus_name = ",".join(genus_names)
        else:
            genus_name = line[1]
        final_dict[qgenome]["genus_names"] = genus_name

        # Species
        final_dict[qgenome]["species_cluster"] = line[2]
        species_names = list(set(spid_ref[line[2]]))
        if species_names:
            if len(species_names) == 1:
                species_name = species_names[0]
            else:
                species_name = ",".join(species_names)
        else:
            species_name = line[2]
        final_dict[qgenome]["species_names"] = species_name


        # get the anicalc results for the genome, sort them by qcov
        genome_df = anicalc_df.loc[qgenome,]
        if len(genome_df.shape) == 1:
            genome_df = genome_df.to_frame().T

        genome_df = genome_df.sort_values(by="qcov", ascending=False)
        genome_df.set_index("tname", inplace=True)

        # iterate the target genomes and store the most similar genome and reference
        final_dict[qgenome]["most_similar_genome"] = genome_df.index[0]
        for tgenome in genome_df.index:
            if tgenome in crass_taxonomy:
                final_dict[qgenome]["most_similar_ref_genome"] = tgenome
                final_dict[qgenome]["ref_genus"] = crass_taxonomy[tgenome]["genus"]
                final_dict[qgenome]["ref_species"] = crass_taxonomy[tgenome]["species"]
                final_dict[qgenome]["pid"] = genome_df.loc[tgenome, "pid"]
                final_dict[qgenome]["qcov"] = genome_df.loc[tgenome, "qcov"]
                break


# ---------------------------------------------------------
# format final_dict as a DataFrame and write to output file
to_write = list()
for genome, columns in final_dict.items():
    to_add = [genome] + [final_dict[genome][column] for column in columns]
    to_write.append(to_add)

to_write_df = pd.DataFrame(to_write, columns=["genome"] + list(columns.keys()))
to_write_df.set_index("genome", inplace=True)
to_write_df.to_csv(snakemake.output[0], index=True, sep="\t")
