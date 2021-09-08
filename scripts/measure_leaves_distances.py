from ete3 import Tree
import pandas as pd

# env: phylogenies.yaml

# read taxonomy file and save it in a dict
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
for line in lines:
    crass_taxonomy[line[0]] = {
                               "family":line[2],
                               "subfamily":line[3],
                               "genus":line[4]
                              }

# read marker tree
t = Tree(snakemake.input[0], format=1)


# assign taxonomy
crAssUS_genomes = list()
for leaf in t.iter_leaves():
    # check if the leaf comes from the reference set
    genome = leaf.name.split("|")[0]
    if genome in crass_taxonomy:
        leaf.add_features(family=crass_taxonomy[genome]["family"],
                          subfamily=crass_taxonomy[genome]["subfamily"],
                          genus=crass_taxonomy[genome]["genus"],
                          genome=genome)
    else:
        leaf.add_features(family="",
                          subfamily="",
                          genus="",
                          genome=genome)
        crAssUS_genomes.append(genome)


# set the dict to store the distances
all_genomes = list(list(crass_taxonomy.keys()) + crAssUS_genomes)
distances = {crAssUS_genome:{genome:"NA" for genome in all_genomes} for crAssUS_genome in crAssUS_genomes}


# traverse the leaves looking for the crAssUS genomes. Once found, measure the distance to the rest of the leaves
for query_leaf in t.iter_leaves():
    if query_leaf.genome in crAssUS_genomes:
        # measure the distance of this branch to the rest of the tree
        for leaf in t.iter_leaves():
            dist = t.get_distance(query_leaf, leaf)
            distances[query_leaf.genome][leaf.genome] = dist


# convert to dataframe and save the matrix to file
df = pd.DataFrame(distances)
df.to_csv(snakemake.output[0], sep="\t")
