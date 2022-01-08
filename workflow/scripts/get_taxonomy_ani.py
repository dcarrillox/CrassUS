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

# get which genids contain more than one ref genus
multiple_ref_gen = list()
for genid, ref_genera in genid_ref.items():
    if ref_genera:
        uniq = list(set(ref_genera))
        if len(uniq) > 1:
            multiple_ref_gen.append(genid)
        else:
            genid_ref[genid] = uniq[0]
    else:
        genid_ref[genid] = genid


# do the same for the species rank
spid_ref = = {line[2]:list() for line in lines}
for line in lines:
    if line[0] in crass_taxonomy:
        spid_ref[line[2]] += [crass_taxonomy[line[0]][1]]

# get which genids contain more than one ref species
multiple_ref_sp = list()
for genid, ref_sps in spid_ref.items():
    if ref_sps:
        uniq = list(set(ref_sps))
        if len(uniq) > 1:
            multiple_ref_sp.append(genid)
        else:
            spid_ref[genid] = uniq[0]
    else:
        spid_ref[genid] = genid


# read anicalc table so I can get later the most similar ref genome
anicalc_df = pd.read_csv(snakemake.input.anicalc[0], sep="\t", header=0, index_col=0)
anicalc_genomes = list(set(anicalc_df.index.tolist()))


# iterate the assignments table again and the anicalc table to get the most similar ref genome, store in a dict
most_similars = {line[0]:{"most_similar_genome":str(), "most_similar_ref":str(), "AF":0} for line in lines}

for line in lines:
    qgenome = line[0]
    # skip reference genomes
    if qgenome not in crass_taxonomy and qgenome in anicalc_genomes:
        # get the anicalc results for the genome, sort them by qcov
        genome_df = anicalc_df.loc[qgenome,]
        if len(genome_df.shape) == 1:
            genome_df = genome_df.to_frame().T

        genome_df = genome_df.sort_values(by="qcov", ascending=False)
        genome_df.set_index("tname", inplace=True)

        # iterate the target genomes and store the most similar genome and reference
        most_similars[qgenome]["most_similar_genome"] = genome_df.index[0]
        for tgenome in genome_df.index:
            if tgenome in crass_taxonomy:
                most_similars[qgenome]["most_similar_ref"] = crass_taxonomy[tgenome][0]
                most_similars[qgenome]["AF"] = genome_df.loc[tgenome, "qcov"]
                break




# iterate the genomes again, this time noticing those genomes within a genid with several ref genera or species
to_write = list()
for line in lines:
    qgenome = line[0]
    if qgenome not in crass_taxonomy:
        to_add = [qgenome]

        if line[1] not in multiple_ref_gen:
            to_add.append(genid_ref[line[1]])
        # if so, add the most similar ref genus
        else:
            print(qgenome, ", multiple ref genus in the genid")
            to_add.append(most_similars[qgenome]["most_similar_ref"])

        to_add.append(most_similars[qgenome]["most_similar_ref"])
        to_add.append(most_similars[qgenome]["AF"])

        # assign species
        if line[2] in spid_ref:
            to_add.append(spid_ref[line[2]])
        else:
            to_add.append(line[2])


        to_write.append(to_add)


to_write_df = pd.DataFrame(to_write, columns=["genome", "genus", "most_similar_ref", "AF", "species"])
to_write_df.to_csv(snakemake.output[0], index=False, sep="\t")
