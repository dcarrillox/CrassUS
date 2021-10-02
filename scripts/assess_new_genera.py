from ete3 import Tree
import pandas as pd
import os

# env: phylogenies.yaml

# read taxonomy assigned by marker trees
markers_taxa = pd.read_csv(snakemake.input.taxa_markers[0], sep="\t", header=0, index_col=0)

# read shared_taxonomy_table and grab the new_genus ones
taxa_shared_table = pd.read_csv(snakemake.input.taxa_shared[0], sep="\t", header=0, index_col=0, keep_default_na=False)
new_genera_table = taxa_shared_table[taxa_shared_table["genera"] == "new_genus"]

# convert to dict
genomes_similars = dict()
for genome, similars in zip(new_genera_table.index, new_genera_table["most_similar_genomes"].astype(str)):
    if similars != "":
        genomes_similars[genome] = similars.split(",")
    else:
        genomes_similars[genome] = []


# create the clusters of new_genus
# init a dict with 'new_genus_0'. It will be removed later
cont = 0
new_genera = {"new_genus_0":list()}
for genome, similars in genomes_similars.items():
    # put genome and similars together in a list
    all_queries = similars + [genome]
    # iterate the genera already in the new_genera dict. If there are no matches,
    # create a new genus
    check = False
    for target_genus, genomes in new_genera.items():
        for query_genome in all_queries:
            # match
            if query_genome in genomes:
                new_genera[target_genus] += all_queries
                check = True
    # no matches
    if not check:
        cont += 1
        new_genera[f"new_genus_{cont}"] = all_queries


# create another dict from the previous one, k=genome  v=new_genus. It will easier
# later when annotating the tree
del new_genera["new_genus_0"]
genomes_new_genus = dict()
for new_genus, genomes in new_genera.items():
    for genome in genomes:
        genomes_new_genus[genome] = list()

for new_genus, genomes in new_genera.items():
    for genome in list(set(genomes)):
        genomes_new_genus[genome] += [new_genus]

# check there is only one genus assigned to each genome
for genome, new_genus in genomes_new_genus.items():
    if len(new_genus) > 1:
        print(f"Watch out! Genome {genome} has been assigned to two different new genera {new_genus}.")


# put the new genera in the trees and look for monophyletic clades
for marker_tree in snakemake.input.markers_trees:
    # get the marker
    marker = os.path.basename(marker_tree).split("_trimmed.nwk")[0]

    # read crass_reference taxonomic classification
    crass_taxonomy = dict()
    lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
    for line in lines:
        crass_taxonomy[line[0]] = {
                                   "family":line[2],
                                   "subfamily":line[3],
                                   "genus":line[4]
                                  }

    # read tree
    t = Tree(marker_tree, format=1)
    # assign taxonomy for the reference crAssphages
    for leaf in t.iter_leaves():
        # check if the leaf comes from the reference set
        genome = leaf.name.split("|")[0]
        if genome in crass_taxonomy:
            leaf.add_features(family=crass_taxonomy[genome]["family"],
                              subfamily=crass_taxonomy[genome]["subfamily"],
                              genus=crass_taxonomy[genome]["genus"],
                              genome=genome)

        # it is not from the reference set
        else:
            # check if the genome is in a new_genus. In that case, set "new_genus_X"
            # as genus. Otherwise keep the genus from the markers table
            if genome in genomes_new_genus:
                leaf.add_features(genus=genomes_new_genus[genome][0], genome=genome)
                # set the new genomes in the markers table, last column (shared)
                markers_taxa.loc[genome, "shared_content"] = genomes_new_genus[genome][0]
            else:
                leaf.add_features(genus=markers_taxa.loc[genome, f"genus_{marker}"], genome=genome)


    # find the LCA of the two outgroup species
    outgs_leaves = t.search_nodes(family="outgroup")
    outgs_lca = t.get_common_ancestor(outgs_leaves)
    # reroot the tree
    t.set_outgroup(outgs_lca)

    # iterate the new_genera looking for monophyletic clades
    for new_genus in new_genera:
        lcas = t.get_monophyletic(values=[new_genus, "unknown"], target_attr="genus")
        # if monophyletic clades containing new_genus and new genomes were found, iterate them
        if lcas:
            # iterate the lcas
            for lca in lcas:
                new_genus_leaves = lca.search_nodes(genus=new_genus)
                # to get the common ancestor, it does not work if there is only one leaf
                if len(new_genus_leaves) > 1:
                    final_lca = lca.get_common_ancestor(new_genus_leaves)
                    # iterate the final_lca while assigning taxonomy
                    for leaf in final_lca.iter_leaves():
                        if leaf.genus == "unknown":
                            #print("\t\tTHIS!!:", leaf.genome, new_genus)
                            markers_taxa.loc[leaf.genome, f"genus_{marker}"] = new_genus
                            # # recalculate completeness, this time based on genus
                            # genome_length = int(leaf.genome.split("_")[-1])
                            # completenes = round((genome_length*100)/taxas_lengths[genus], 2)
                            # markers_taxa.loc[leaf.genome, "completeness"] = completenes
                            # markers_taxa.loc[leaf.genome, "reference"] = genus

markers_taxa.to_csv(snakemake.output[0], index=True, sep="\t")
