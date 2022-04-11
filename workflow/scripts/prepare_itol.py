#!/usr/bin/env python
# env: phylogenies.yaml

from ete3 import Tree
import os


os.makedirs(snakemake.params.dir, exist_ok=True)

families_colors = {"Intestiviridae":"#EE3B3B",
                      "Crevaviridae":"#EE9A00",
                      "Suoliviridae":"#4169E1",
                      "Steigviridae":"#00CED1",
                      "Epsilon":"#CD2990",
                      "Zeta":"#006400",
                      "outgroup": "gray",
                      "NA":"gray"}

outgroups = {"TerL":"NC_021803|775|90",
             "MCP":"NC_021803|443|97",
             "portal":"NC_021803|812|91"}



# ------------------------
# Parse reference taxonomy
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()[1:]]
for line in lines:
    crass_taxonomy[line[0]] = {"family":line[1],
                               "subfamily":line[2],
                               "genus":line[3],
                               "species":line[4]}

# get marker and adjust taxonomy for Akihdevirus balticus (KC821624) with TerL
marker = os.path.basename(snakemake.input[0]).split("_trimmed.nwk")[0]
if marker == "TerL":
    crass_taxonomy["KC821624_1_100418"] = {"family":"NA",
                                           "subfamily":"NA",
                                           "genus":"NA",
                                           "species":"Akihdevirus balticus"}
else:
    crass_taxonomy["KC821624_1_100418"] = {"family":"Steigviridae",
                                           "subfamily":"Asinivirinae",
                                           "genus":"Akihdevirus",
                                           "species":"Akihdevirus balticus"}


# ----------
# Parse tree
t = Tree(snakemake.input[0], format=1)
# root with outgroup
out_leaf = t.search_nodes(name=outgroups[marker])[0]
t.set_outgroup(out_leaf)
# iterate branches. For reference genomes, replace label by its species name and
# save its iToL annotation in a list
leaves_annotation = list()
species_done = set()
for leaf in t.iter_leaves():
    # check if the leaf comes from the reference set
    genome = leaf.name.split("|")[0]
    if genome in crass_taxonomy:
        leaf.add_features(family=crass_taxonomy[genome]["family"])
        species = crass_taxonomy[genome]["species"].replace(" ", "_")
        # check if the name has been already used. If so, try adding "_" until
        # it becomes a new species id
        if f"{species}___1" in species_done:
            for i in range(2,1000):
                tmp = f"{species}___{i}"
                if tmp not in species_done:
                    species_done.add(tmp)
                    species = tmp
                    break
        else:
            species = f"{species}___1"
            species_done.add(species)

        leaf.name = species
        # create annot line for the leaf with the following structure:
        #ID,TYPE,WHAT,COLOR,WIDTH_OR_SIZE_FACTOR,STYLE,BACKGROUND_COLOR
        color = families_colors[crass_taxonomy[genome]["family"]]
        leaf_annot = f"{species}\tlabel\tnode\t{color}\t2\tbold-italic"
        leaves_annotation.append(leaf_annot)

# find LCA of the families and color all their child branches
for family in families_colors.keys():
    if family not in ["outgroup", "NA"]:
        family_leaves = t.search_nodes(family=family)
        family_lca = t.get_common_ancestor(family_leaves)
        # get name and create its annotation
        family_lca.name = family
        color = families_colors[family]
        lca_annot = f"{family}\tbranch\tclade\t{color}"
        leaves_annotation.append(lca_annot)



# write output rooted tree with reference genome labels changed to species name
t.write(format=1, outfile=snakemake.output.tree)


# ----------------------------
# Prepare iToL annotation file
legend_colors = '\t'.join(list(families_colors.values()))
legend_families = '\t'.join(list(families_colors.keys()))
itol_annotation = ["DATASET_STYLE",
                   "SEPARATOR\tTAB",
                   f"DATASET_LABEL\tFamilies {marker}",
                   "COLOR\t#ffff00",
                   f"LEGEND_TITLE\tFamilies {marker}",
                   f"LEGEND_COLORS\t{legend_colors}",
                   f"LEGEND_LABELS\t{legend_families}",
                   "LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\t1\t1",
                   "DATA"]
# include leaves annotation
final_annot_file = itol_annotation + leaves_annotation

with open(snakemake.output.annot, "w") as fout:
    for line in final_annot_file:
        fout.write(line + "\n")
