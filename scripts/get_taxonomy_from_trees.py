from ete3 import Tree
import pandas as pd
import os

# env: phylogenies.yaml

# list the contigs that were found by crAssUS.
df = pd.read_csv(snakemake.input.markers_summary, sep="\t")
crassus_contigs = df["contig"].to_list()

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
        lcas = t.get_monophyletic(values=[subfamily, "new"], target_attr="subfamily")
        # if monophyletic clades containing subfamily and new genomes were found, iterate them
        if lcas:
            # iterate the lcas
            for lca in lcas:
                subfamily_leaves = lca.search_nodes(subfamily=subfamily)
                # to get the common ancestor, it does not work if there is only one leaf
                if len(subfamily_leaves) > 1:
                    final_lca = lca.get_common_ancestor(subfamily_leaves)
                    # iterate the final_lca while assigning taxonomy
                    for leaf in final_lca.iter_leaves():
                        if leaf.subfamily == "new":
                            crassus_classification[leaf.genome][f"subfamily_{marker}"] = subfamily

    for genus in genera:
        lcas = t.get_monophyletic(values=[genus, "new"], target_attr="genus")
        # if monophyletic clades containing genus and new genomes were found, iterate them
        if lcas:
            # iterate the lcas
            for lca in lcas:
                genus_leaves = lca.search_nodes(genus=genus)
                # to get the common ancestor, it does not work if there is only one leaf
                if len(genus_leaves) > 1:
                    final_lca = lca.get_common_ancestor(genus_leaves)
                    # iterate the final_lca while assigning taxonomy
                    for leaf in final_lca.iter_leaves():
                        if leaf.genus == "new":
                            crassus_classification[leaf.genome][f"genus_{marker}"] = genus





    # for subfam in subfamilies:
    # # get all the leaves
    #     subfam_leaves = t.search_nodes(subfamily=subfam)
    #     if len(subfam_leaves) > 1:
    #         subfam_lca = t.get_common_ancestor(subfam_leaves)
    #         # check if it is a monophyletic clade
    #         if not check_monophyletic_subfamily(subfam_lca, subfam):
    #             # name inner nodes
    #             cont = 0
    #             for node in subfam_lca.traverse():
    #                 if not node.is_leaf():
    #                     node.name = f"node_{cont}"
    #                     cont += 1
    #             # get the monophyletic clades
    #             monophyletic_clades = get_monophyletic_subfamily(subfam_lca, subfam)
    #         else:
    #             monophyletic_clades = subfam_lca
    #
    #         # iterate the monophyletic clades while assigning taxonomy
    #         for monophyletic_clade in monophyletic_clades:
    #             for leaf in monophyletic_clade.iter_leaves():
    #                 if leaf.subfamily == "new":
    #                     crassus_classification[leaf.genome][f"subfamily_{marker}"] = subfam
    #
    # for genus in genera:
    # # get all the leaves
    #     genus_leaves = t.search_nodes(genus=genus)
    #     if len(genus_leaves) > 1:
    #         genus_lca = t.get_common_ancestor(genus_leaves)
    #         # check if it is a monophyletic clade
    #         if not check_monophyletic_genus(genus_lca, genus):
    #             # name inner nodes
    #             cont = 0
    #             for node in genus_lca.traverse():
    #                 if not node.is_leaf():
    #                     node.name = f"node_{cont}"
    #                     cont += 1
    #             # get the monophyletic clades
    #             monophyletic_clades = get_monophyletic_genus(genus_lca, genus)
    #         else:
    #             monophyletic_clades = genus_lca
    #
    #         # iterate the monophyletic clades while assigning taxonomy
    #         for monophyletic_clade in monophyletic_clades:
    #             for leaf in monophyletic_clade.iter_leaves():
    #                 if leaf.genus == "new":
    #                     crassus_classification[leaf.genome][f"genus_{marker}"] = genus


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
