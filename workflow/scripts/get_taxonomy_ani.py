import pandas as pd

# envs: compare_genomes.yaml


# read crass taxonomy, store genus and species
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()[1:]]
crass_taxonomy = dict()
for line in lines:
    crass_taxonomy[line[0]] = [line[3], line[4]]


# read ani_clusters file and store the genus__XX that correspond to a reference genus
lines = [line.strip().split("\t") for line in open(snakemake.input.assignments[0]).readlines()[1:]]
genid_ref = {line[1]:list() for line in lines}
for line in lines:
    if line[0] in crass_taxonomy:
        genid_ref[line[1]] += [crass_taxonomy[line[0]][0]]

# do the same for the species rank
spid_ref = {line[2]:list() for line in lines}
for line in lines:
    if line[0] in crass_taxonomy:
        spid_ref[line[2]] += [crass_taxonomy[line[0]][1]]

# read anicalc table so I can get later the most similar ref genome
anicalc_df = pd.read_csv(snakemake.input.anicalc[0], sep="\t", header=0, index_col=0)
anicalc_genomes = set(anicalc_df.index.tolist())


# iterate the assignments table again and the anicalc table to get the most similar ref genome, store in a dict
final_dict = {line[0]:
                        {"genus_cluster":str(),
                        "genus_names":str(),
                        "species_cluster":str(),
                        "species_names":str(),
                        "most_similar_genome":str(),
                        "most_similar_ref":str(),
                        "ref_genus":str(),
                        "ref_species":str(),
                        "pid":str(),
                        "qcov":0}
                for line in lines if line[0] not in crass_taxonomy and line[0] in anicalc_genomes}

for line in lines:
    qgenome = line[0]

    if line[0] in final_dict:

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
                final_dict[qgenome]["most_similar_ref"] = tgenome
                final_dict[qgenome]["ref_genus"] = crass_taxonomy[tgenome][0]
                final_dict[qgenome]["ref_species"] = crass_taxonomy[tgenome][1]
                final_dict[qgenome]["pid"] = genome_df.loc[tgenome, "pid"]
                final_dict[qgenome]["qcov"] = genome_df.loc[tgenome, "qcov"]
                break





# # iterate the genomes again, this time noticing those genomes within a genid with several ref genera or species
# to_write = list()
# for line in lines:
#     qgenome = line[0]
#     if qgenome not in crass_taxonomy:
#         to_add = [qgenome]
#
#         if line[1] not in multiple_ref_gen:
#             to_add.append(genid_ref[line[1]])
#         # if so, add the most similar ref genus
#         else:
#             print(qgenome, ", multiple ref genus in the genid")
#             to_add.append(most_similars[qgenome]["most_similar_ref"])
#
#         to_add.append(most_similars[qgenome]["most_similar_ref"])
#         to_add.append(most_similars[qgenome]["AF"])
#
#         # assign species
#         if line[2] in spid_ref:
#             to_add.append(spid_ref[line[2]])
#         else:
#             to_add.append(line[2])
#
#
#         to_write.append(to_add)


#to_write_df = pd.DataFrame(to_write, columns=["genome", "genus", "most_similar_ref", "AF", "species"])
#to_write_df.to_csv(snakemake.output[0], index=False, sep="\t")

to_write = list()
for genome, columns in final_dict.items():
    to_add = [genome] + [final_dict[genome][column] for column in columns]
    to_write.append(to_add)

to_write_df = pd.DataFrame(to_write, columns=["genome"] + list(columns.keys()))
to_write_df.set_index("genome", inplace=True)
to_write_df.to_csv(snakemake.output[0], index=True, sep="\t")
