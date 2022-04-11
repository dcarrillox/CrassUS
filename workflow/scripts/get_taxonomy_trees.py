from ete3 import Tree
import pandas as pd
import os

# env: phylogenies.yaml

# --------------------------------------------
# list the contigs that were found by crAssUS.
markers_summary = pd.read_csv(snakemake.input.markers_summary, sep="\t", header=0, index_col=0)
crassus_contigs = markers_summary.index.to_list()



# ----------------------------------------
# get which markers underwent the analysis
markers = [os.path.basename(tree_file).split("_trimmed.nwk")[0] for tree_file in snakemake.input.markers_trees if "iToL" not in tree_file]
markers = sorted(markers, reverse=True)



# ------------------------------------------------------------------------------
# init the taxa classification dict, containing all the markers that underwent the analysis
# init them as unknown, if they were classified (or not) the value will be replaced
crassus_classification = {contig:dict() for contig in crassus_contigs}
for contig in crassus_classification:
    for rank in ["family", "subfamily", "genus"]:
        for marker in markers:
            crassus_classification[contig][f"{rank}_{marker}"] = "unknown"
            crassus_classification[contig][f"{rank}_{marker}"] = "unknown"
            crassus_classification[contig][f"{rank}_{marker}"] = "unknown"



# ---------------------------------------------
# read crass_reference taxonomic classification
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()[1:]]
for line in lines:
    crass_taxonomy[line[0]] = {"family":line[1],"subfamily":line[2],"genus":line[3]}

# get all the possible names for the three ranks
families = list(set([crass_taxonomy[genome]["family"] for genome in crass_taxonomy]))
subfamilies = list(set([crass_taxonomy[genome]["subfamily"] for genome in crass_taxonomy]))
genera = list(set([crass_taxonomy[genome]["genus"] for genome in crass_taxonomy]))

# remove outgroup and NA
families.remove("outgroup")
subfamilies.remove("outgroup")
subfamilies.remove("")
genera.remove("outgroup")




# -------------------------------------------
# go throug the markers trees and parse them
outgroups = {"TerL":"NC_021803|775|90", "MCP":"NC_021803|443|97", "portal":"NC_021803|812|91"}

# discard iToL files from the input list
markers_trees =[file for file in snakemake.input.markers_trees if "iToL" not in file]

for marker_tree in markers_trees:
    # get the marker
    marker = os.path.basename(marker_tree).split("_trimmed.nwk")[0]


    # correct Akihdevirus balticus (KC821624) for TerL
    if marker == "TerL":
        crass_taxonomy["KC821624_1_100418"] = {"family":"NA",
                                               "subfamily":"NA",
                                               "genus":"NA"}
    else:
        crass_taxonomy["KC821624_1_100418"] = {"family":"Steigviridae",
                                               "subfamily":"Asinivirinae",
                                               "genus":"Akihdevirus"}


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

    # find the LCA of the outgroup species and reroot with it
    out_leaf = t.search_nodes(name=outgroups[marker])[0]
    t.set_outgroup(out_leaf)


    # -------------------------
    # classify crassus_contigs

    # Family rank
    for family in families:
        family_leaves = t.search_nodes(family=family)
        family_lca = t.get_common_ancestor(family_leaves)
        # iterate the leaves of the LCA
        for leaf in family_lca.iter_leaves():
            if leaf.family == "new":
                crassus_classification[leaf.genome][f"family_{marker}"] = family


    # Subfamily rank
    for subfamily in subfamilies:
        lcas = t.get_monophyletic(values=[subfamily, "new"], target_attr="subfamily")
        # if monophyletic clades containing subfamily and new genomes were found, iterate them
        if lcas:
            for lca in lcas:
                subfamily_leaves = lca.search_nodes(subfamily=subfamily)
                # to get the common ancestor, it doesn't work if there is only one leaf
                if len(subfamily_leaves) > 1:
                    final_lca = lca.get_common_ancestor(subfamily_leaves)
                    # iterate the final_lca while assigning taxonomy
                    for leaf in final_lca.iter_leaves():
                        if leaf.subfamily == "new":
                            crassus_classification[leaf.genome][f"subfamily_{marker}"] = subfamily


    # Genus rank
    for genus in genera:
        lcas = t.get_monophyletic(values=[genus, "new"], target_attr="genus")
        # if monophyletic clades containing genus and new genomes were found, iterate them
        if lcas:
            for lca in lcas:
                genus_leaves = lca.search_nodes(genus=genus)
                # to get the common ancestor, it doesn't work if there is only one leaf
                if len(genus_leaves) > 1:
                    final_lca = lca.get_common_ancestor(genus_leaves)
                    # iterate the final_lca while assigning taxonomy
                    for leaf in final_lca.iter_leaves():
                        if leaf.genus == "new":
                            crassus_classification[leaf.genome][f"genus_{marker}"] = genus

    # check which genomes are not present in the tree and write the reason
    # according to "markers_summary.txt"
    for genome in crassus_contigs:
        if genome not in genomes_marker:
            # get what happened in this genome from the summary_markers file
            reason = markers_summary.loc[genome, marker]
            crassus_classification[genome][f"family_{marker}"] = reason
            crassus_classification[genome][f"subfamily_{marker}"] = reason
            crassus_classification[genome][f"genus_{marker}"] = reason



# ---------------------------------------------------------------
# convert the classification dictionary to table and write to csv
df = pd.DataFrame(crassus_classification)
tdf = df.transpose()
tdf.to_csv(snakemake.output[0], sep="\t", index=True)
