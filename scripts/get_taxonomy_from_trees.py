from ete3 import Tree
import pandas as pd
import os

# list the contigs that were found by crAssUS.
# crassus_contigs = [os.path.basename(file).replace(".fasta", "") for file in snakemake.input.crassus_contigs]
# print(crassus_contigs)

#print(snakemake.input.markers_summary)
df = pd.read_csv(snakemake.input.markers_summary, sep="\t")
crassus_contigs = df["contig"].to_list()
#print(crassus_contigs)

# get which markers underwent the analysis
print(snakemake.input.markers_trees)
markers = [os.path.basename(tree_file).split("_trimmed.nwk")[0] for tree_file in snakemake.input.markers_trees]
markers = sorted(markers, reverse=True)

# init the taxa classification dict, containing all the markers that underwent the analysis
crassus_classification = {contig:dict() for contig in crassus_contigs}
for contig in crassus_classification:
    for marker in markers:
        crassus_classification[contig][f"family_{marker}"] = "unknown"
        crassus_classification[contig][f"subfamily_{marker}"] = "unknown"
        crassus_classification[contig][f"genus_{marker}"] = "unknown"


# go through the markers trees and parse them
for marker_tree in snakemake.input.markers_trees:
    # get the marker
    marker = os.path.basename(marker_tree).split("_trimmed.nwk")[0]
    print(marker)

    # read crass_reference taxonomic classification
    crass_taxonomy = dict()
    lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
    for line in lines:
        crass_taxonomy[line[0]] = {
                                   "family":line[2],
                                   "subfamily":line[3],
                                   "genus":line[4]
                                  }

    # read the file with the monophyletic fix for the taxonomy
    fix_file = f"resources/{marker}_fixed_monophyl_tax.txt"
    print(fix_file)
    lines = [line.strip().split("\t") for line in open(fix_file).readlines()[1:]]  # discard header
    # go through the lines and store the change in the taxa dict from above
    for line in lines:
        genome, new_subfam, new_gen = line[0], line[2], line[4]

        if new_subfam != "":
            crass_taxonomy[genome]["subfamily"] = new_subfam
        if new_gen != "":
            crass_taxonomy[genome]["genus"] = new_gen

    # init a list to store which genomes were included in the tree, so I can
    # say "Not found" in the final table
    genomes_marker = list()

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
        else:
            genomes_marker.append(genome)
            leaf.add_features(family="new",
                              subfamily="new",
                              genus="new",
                              genome=genome)

    # find the LCA of the two outgroup species
    outgs_leaves = t.search_nodes(family="outgroup")
    outgs_lca = t.get_common_ancestor(outgs_leaves)
    # reroot the tree
    t.set_outgroup(outgs_lca)


    ## classify crassus_contigs
    # first get all the possible names for the three ranks
    families = list(set([crass_taxonomy[genome]["family"] for genome in crass_taxonomy]))
    subfamilies = list(set([crass_taxonomy[genome]["subfamily"] for genome in crass_taxonomy]))
    genera = list(set([crass_taxonomy[genome]["genus"] for genome in crass_taxonomy]))


    for family in families:
        family_leaves = t.search_nodes(family=family)
        family_lca = t.get_common_ancestor(family_leaves)
        # iterate the leaves of the LCA
        for leaf in family_lca.iter_leaves():
            if leaf.family == "new":
                crassus_classification[leaf.genome][f"family_{marker}"] = family

    for subfamily in subfamilies:
        subfamily_leaves = t.search_nodes(subfamily=subfamily)
        if len(subfamily_leaves) > 1:
            subfamily_lca = t.get_common_ancestor(subfamily_leaves)
            # iterate the leaves of the LCA
            for leaf in subfamily_lca.iter_leaves():
                if leaf.subfamily == "new":
                    crassus_classification[leaf.genome][f"subfamily_{marker}"] = subfamily

    for genus in genera:
        genus_leaves = t.search_nodes(genus=genus)
        if len(genus_leaves) > 1:
            genus_lca = t.get_common_ancestor(genus_leaves)
            # iterate the leaves of the LCA
            for leaf in genus_lca.iter_leaves():
                if leaf.genus == "new":
                    crassus_classification[leaf.genome][f"genus_{marker}"] = genus




    # # classify crassus contigs
    # for leaf in t.iter_leaves():
    #     # check if it is a new contig found by crAssUS
    #     if leaf.genome in crassus_classification:
    #
    #         # family and subfamily: inspect the upper node
    #         node = leaf.up
    #         fam,subfam = list(), list()
    #         for hoja in node.iter_leaves():
    #             fam.append(hoja.family)
    #             subfam.append(hoja.subfamily)
    #
    #         fam.remove("")
    #         subfam.remove("")
    #
    #         fam, subfam = list(set(fam)), list(set(subfam))
    #         if len(fam) == 1:
    #             crassus_classification[leaf.genome][f"family_{marker}"] = fam[0]
    #         if len(subfam) == 1:
    #             crassus_classification[leaf.genome][f"subfamily_{marker}"] = subfam[0]
    #
    #         # genus: inspect the two upper nodes
    #         node = leaf.up.up
    #         genus = list()
    #         for hoja in node.iter_leaves():
    #             genus.append(hoja.genus)
    #
    #         genus.remove("")
    #         genus = list(set(genus))
    #         if len(genus) == 1:
    #             crassus_classification[leaf.genome][f"genus_{marker}"] = genus[0]


    # check which genomes are not present in the tree and call them as "Not found"
    # for that marker
    for genome in crassus_contigs:
        if genome not in genomes_marker:
            crassus_classification[genome][f"family_{marker}"] = "Not found"
            crassus_classification[genome][f"subfamily_{marker}"] = "Not found"
            crassus_classification[genome][f"genus_{marker}"] = "Not found"


# convert the dictionary to table and write to csv with pandas
#print(crassus_classification)
df = pd.DataFrame(crassus_classification)
tdf = df.transpose()
tdf.to_csv(snakemake.output[0], sep="\t", index=True)
