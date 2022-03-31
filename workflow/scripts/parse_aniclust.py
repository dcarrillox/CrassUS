
# read reference taxonomy
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
for line in lines:
    crass_taxonomy[line[0]] = {
                               "family":line[2],
                               "subfamily":line[3],
                               "genus":line[4]
                              }


# read aniclust genus file. Init a dict to store the putative genus for all the
# genomes that are not in the reference, ie. genomes found by crassus
crassus_genomes_genus = dict()
# retain the column [1] with all the genomes and split by ","
clusters = [line.strip().split("\t")[1].split(",") for line in open(snakemake.input.genus).readlines()]
cont = 1
for cluster in clusters:
    # init the ref genome string
    ref_genera = list()
    crassus_genomes = list()
    for genome in cluster:
        if genome in crass_taxonomy:
            ref_genera.append(crass_taxonomy[genome]["genus"])
        else:
            crassus_genomes.append(genome)

    # if ref genomes were present in the cluster
    if ref_genera:
        genus = list(set(ref_genera))
        # check there is only one genus
        if len(genus) > 1:
            print(f"Genus discrepancies {genus}\n\tcluster:{cluster}\n")
        else:
            for genome in crassus_genomes:
                crassus_genomes_genus[genome] = genus[0]

    else:
        genus = f'new_genus_{cont}'
        for genome in crassus_genomes:
            crassus_genomes_genus[genome] = genus
        cont += 1

# write genus assignments
with open(snakemake.output.genus, "w") as fout:
    for genome, genus in crassus_genomes_genus.items():
        fout.write(f"{genome}\t{genus}\n")



# read aniclust genus file. Init a dict to store the putative genus for all the
# genomes that are not in the reference, ie. genomes found by crassus
crassus_genomes_species = dict()
# retain the column [1] with all the genomes and split by ","
clusters = [line.strip().split("\t")[1].split(",") for line in open(snakemake.input.species).readlines()]
cont = 1
for cluster in clusters:
    # init the ref genome string
    ref_species = list()
    crassus_genomes = list()
    for genome in cluster:
        if genome in crass_taxonomy:
            ref_species.append(genome)
        else:
            crassus_genomes.append(genome)

    # if ref genomes were present in the cluster
    if ref_species:
        for genome in crassus_genomes:
            crassus_genomes_species[genome] = "present in ref set"

    else:
        species = f'new_species_{cont}'
        for genome in crassus_genomes:
            crassus_genomes_species[genome] = species
        cont += 1

# write genus assignments
with open(snakemake.output.species, "w") as fout:
    for genome, species in crassus_genomes_species.items():
        fout.write(f"{genome}\t{species}\n")
