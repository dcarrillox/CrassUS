from ete3 import Tree, TreeStyle, TextFace, RectFace, NodeStyle
import os

# list the contigs that were found by crAssUS
crassus_contigs = [os.path.basename(file).replace(".fasta", "") for file in snakemake.input.found_contigs]
# init a dict for their classification
crassus_terl_classification = {contig:{"family_terl":str(), "subfamily_terl":str(), "genus_terl":str()} for contig in crassus_contigs}


# read taxonomic classification for the reference crAssphages
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
for line in lines:
    crass_taxonomy[line[0]] = {
                               "family":line[2],
                               "subfamily":line[3],
                               "genus":line[4]
                              }


# read tree
t = Tree(snakemake.input.terl_tree[0], format=1)
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
        leaf.add_features(family="",
                          subfamily="",
                          genus="",
                          genome=genome)
# find the LCA of the two outgroup species
outgs_leaves = t.search_nodes(family="outgroup")
outgs_lca = t.get_common_ancestor(outgs_leaves)
# reroot the tree
t.set_outgroup(outgs_lca)


# classify crassus contigs
for leaf in t.iter_leaves():
    # check if it is a new contig found by crAssUS
    if leaf.genome in crassus_terl_classification:

        # family and subfamily: inspect the upper node
        node = leaf.up
        fam,subfam = list(), list()
        for hoja in node.iter_leaves():
            fam.append(hoja.family)
            subfam.append(hoja.subfamily)

        fam.remove("")
        subfam.remove("")

        fam, subfam = list(set(fam)), list(set(subfam))
        if len(fam) == 1:
            crassus_terl_classification[leaf.genome]["family_terl"] = fam[0]
        if len(subfam) == 1:
            crassus_terl_classification[leaf.genome]["subfamily_terl"] = subfam[0]

        # genus: inspect the two upper nodes
        node = leaf.up.up
        genus = list()
        for hoja in node.iter_leaves():
            genus.append(hoja.genus)

        genus.remove("")
        genus = list(set(genus))
        if len(genus) == 1:
            crassus_terl_classification[leaf.genome]["genus_terl"] = genus[0]


with open(snakemake.output[0], "w") as fout:
    for genome, taxa in crassus_terl_classification.items():
        fout.write(f'{genome}\t{taxa["family_terl"]}\t{taxa["subfamily_terl"]}\t{taxa["genus_terl"]}\n')
