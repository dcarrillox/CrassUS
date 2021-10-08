import pandas as pd
import os

# env: utils.yaml

# start by reading the table with marker and shared annotation
annot_table = pd.read_csv(snakemake.input.taxa_table[0], sep="\t", header=0, index_col=0)
# retain only complete genomes
annot_table_complete = annot_table[annot_table["completeness"] >= 90]


# consider only crassphages from the reference and complete crasuss to define species
# read crass_reference taxonomic classification
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
for line in lines:
    crass_taxonomy[line[0]] = {
                               "family":line[2],
                               "subfamily":line[3],
                               "genus":line[4]
                              }
# list with reference and complete crassus
species_complete_genomes = list(crass_taxonomy.keys()) + annot_table_complete.index.tolist()


# create a dict for each of the genus in shared_content.
genera_species = {genus:dict() for genus in annot_table_complete["shared_content"]}
genera_species_count = {genus:0 for genus in annot_table_complete["shared_content"]}
# init a dict for later with k=genome  v=species (species to fill in later)
genomes_species = {genome:list() for genome in annot_table_complete.index}


# iterate the complete genomes and check their pyani results
complete_genomes = annot_table_complete.index
for query_genome in complete_genomes:
    # set the paths to the pyani files
    pyani_cov_file = f"results/7_ANI/1_most_similar/{query_genome}/ANIb_alignment_coverage.tab"
    pyani_sim_file = f"results/7_ANI/1_most_similar/{query_genome}/ANIb_percentage_identity.tab"

    # read pyani files
    cov_df = pd.read_csv(pyani_cov_file, sep="\t", header=0, index_col=0)
    sim_df = pd.read_csv(pyani_sim_file, sep="\t", header=0, index_col=0)

    # check the coverage for the query genome
    target_genomes_cov = list()
    for target_genome in cov_df:
        # if the coverage is >= 0.85, add the target genome to the list
        if cov_df.loc[query_genome, target_genome] >= 0.85 and query_genome != target_genome and target_genome in species_complete_genomes:
            target_genomes_cov.append(target_genome)
            # check if the coverage match is reciprocal
            if not cov_df.loc[target_genome, query_genome] >= 0.85:
                pass
                #print(f"Watch out!! Genome {query_genome} does not cover 85% of {target_genome}.")

    # check the similarity with the query genome
    target_genomes_sim = list()
    for target_genome in sim_df:
        # if the coverage is >= 0.85, add the target genome to the list
        if sim_df.loc[query_genome, target_genome] >= 0.95 and query_genome != target_genome and target_genome in species_complete_genomes:
            target_genomes_sim.append(target_genome)
            # check if the similarity match is reciprocal
            if not sim_df.loc[target_genome, query_genome] >= 0.95:
                pass
                #print(f"Watch out!! Genome {query_genome} is not 95% similar to {target_genome}.")

    # find the intersection between the cov and simi lists. First check that both list are non empty
    cov_sim_target_genomes = list()
    if target_genomes_cov and target_genomes_sim:
        cov_sim_target_genomes = [genome for genome in target_genomes_cov if genome in target_genomes_sim]

    # there were genomes passing cov and sim cutoffs
    if cov_sim_target_genomes:
        # check that all the genomes are annotated with the same genus
        genomes_genus = [annot_table_complete.loc[query_genome, "shared_content"]]
        for genome in cov_sim_target_genomes:
            if genome in crass_taxonomy:
                genomes_genus.append(crass_taxonomy[genome]["genus"])
            else:
                genomes_genus.append(annot_table_complete.loc[query_genome, "shared_content"])
        genomes_genus = list(set(genomes_genus))

        # check there is only one genus associated
        if len(genomes_genus) > 1:
            print(f"Different genera {genomes_genus} matching to {query_genome}. Skipping")
        else:
            genomes_genus = genomes_genus[0]
            # put query and target genomes in the same list
            all_genomes = cov_sim_target_genomes + [query_genome]


            # iterate the different keys of the species dict
            check = True
            for species in genera_species[genomes_genus]:
                # iterate the genomes to check in any is already in a species bin
                for genome in all_genomes:
                    if genome in list(set(genera_species[genomes_genus][species])):
                        genera_species[genomes_genus][species] += all_genomes
                        check = False
                        if "Piconjevirus" in species:
                            print(query_genome, genome, species, genera_species[genomes_genus][species])


            # if any genome was found in a species bin
            if check:
                genera_species_count[genomes_genus] += 1
                cont = genera_species_count[genomes_genus]
                genera_species[genomes_genus][f"{genomes_genus}_sp{cont}"] = all_genomes
                if "Piconjevirus" in genomes_genus:
                    print(query_genome, f"{genomes_genus}_sp{cont}", genera_species[genomes_genus][f"{genomes_genus}_sp{cont}"])
    else:
        ## add single genomes
        genus = annot_table_complete.loc[query_genome, "shared_content"]
        genera_species_count[genus] += 1
        cont = genera_species_count[genus]
        genera_species[genus][f"{genus}_sp{cont}"] = [query_genome]




# convert the dictionary to the form k=genome  v=species
genomes_species = {genome:list() for genome in annot_table_complete.index}
for genus in genera_species:
    for species in genera_species[genus]:
        if "Piconjevirus" in species:
            print(species, list(set(genera_species[genus][species])))
        for genome in list(set(genera_species[genus][species])):
            if genome in genomes_species:
                genomes_species[genome] += [species]



# check only one species was associated with the genomes. If so, remove it from the
# dictionary
to_del = list()
for genome, species in genomes_species.items():
    if len(species) > 1:
        #print(f"More than one species {species} were assigned to {genome}. Skipping.")
        to_del.append(genome)

for genome in to_del:
    del genomes_species[genome]



# put species info into final table
for genome, species in genomes_species.items():
    annot_table.loc[genome, "species"] = species[0]


# write final table
#annot_table.to_csv(snakemake.output[0], index=True, sep="\t")
