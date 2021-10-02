from ete3 import Tree
import pandas as pd
import os

# env: phylogenies.yaml

# read taxa assignments by shared proteins
shared_taxa = pd.read_csv(snakemake.input.taxa_shared[0], sep="\t", header=0, index_col=0)
# keep only complete genomes
shared_taxa = shared_taxa[shared_taxa["complete"] == "yes"]

# read genera's lengths to recalculate later the coverage of the newly classified genomes
lines = [line.strip().split() for line in open(snakemake.params.lengths).readlines()]
taxas_lengths = {line[0]:float(line[1]) for line in lines}


# read taxonomy assigned by marker trees
markers_taxa = pd.read_csv(snakemake.input.taxa_markers[0], sep="\t", header=0, index_col=0)

# iterate the trees
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


    # get the genera
    genera = list(set([crass_taxonomy[genome]["genus"] for genome in crass_taxonomy]))

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
            # get annot from the marker tree
            marker_taxa_genome = markers_taxa.loc[genome, f"genus_{marker}"]

            # check if the genome is complete and was annotated with any other genus
            if genome in shared_taxa.index:
                shared_taxa_genome = shared_taxa.loc[genome, "most_similar_genus"]

                # shared and markers agree
                if shared_taxa_genome == marker_taxa_genome:
                    genus = shared_taxa_genome
                    # add the shared genus to the markers table, it will be the output
                    markers_taxa.loc[genome, "shared_content"] = genus

                # shared new_genus, but genus by marker. In this case probably the
                # shared value is close to 0.7
                elif shared_taxa_genome == "new_genus" and marker_taxa_genome != "unknown":
                    print(f"{genome}: shared value below 0.7, but classified as {marker_taxa_genome} by {marker}.")
                    genus = marker_taxa_genome

                # the other way around, shared genus but not by marker
                # THIS IS THE USEFUL CASE, IT EXPANDS THE MONOPHYLETIC CLADE
                elif shared_taxa_genome != "new_genus" and marker_taxa_genome == "unknown":
                    genus = shared_taxa_genome
                    print("\tLook at:", genome, genus)
                    # add the shared genus to the markers table, it will be the output
                    markers_taxa.loc[genome, "shared_content"] = genus

                # shared and marker genus, but they don't agree. Most likely is that
                # these are two similar genera and the one assigned by marker is a bit
                # below 0.7 shared. Keep marker classification.
                elif shared_taxa_genome != "new_genus" and marker_taxa_genome != "unknown":
                    print(f"{genome}: shared and marker ({marker}) don't agree [{shared_taxa_genome}, {marker_taxa_genome}]. Keeping the marker one, {marker_taxa_genome}.")
                    genus = marker_taxa_genome


                # CHECK: new genus
                elif shared_taxa_genome == "new_genus" and marker_taxa_genome == "unknown":
                    genus = "new_genus"
                    # add the shared genus to the markers table, it will be the output
                    markers_taxa.loc[genome, "shared_content"] = "new_genus"

                else:
                    print(f"\t??? Check this genome: {genome}")

                # add the genus information to the leaf
                leaf.add_features(genus=genus, genome=genome)

            else:
                leaf.add_features(genus=marker_taxa_genome, genome=genome)


    # find the LCA of the two outgroup species
    outgs_leaves = t.search_nodes(family="outgroup")
    outgs_lca = t.get_common_ancestor(outgs_leaves)
    # reroot the tree
    t.set_outgroup(outgs_lca)

    for genus in genera:
        print(genus)
        lcas = t.get_monophyletic(values=[genus, "unknown"], target_attr="genus")
        # if monophyletic clades containing genus and new genomes were found, iterate them
        if lcas:
            print(lcas)
            # iterate the lcas
            for lca in lcas:
                genus_leaves = lca.search_nodes(genus=genus)
                # to get the common ancestor, it does not work if there is only one leaf
                if len(genus_leaves) > 1:
                    final_lca = lca.get_common_ancestor(genus_leaves)
                    # iterate the final_lca while assigning taxonomy
                    for leaf in final_lca.iter_leaves():
                        if leaf.genus == "unknown":
                            #pass
                            print("\t\tTHIS!!:", leaf.genome, genus)
                            markers_taxa.loc[leaf.genome, f"genus_{marker}"] = genus
                            # recalculate completeness, this time based on genus
                            genome_length = int(leaf.genome.split("_")[-1])
                            completenes = round((genome_length*100)/taxas_lengths[genus], 2)
                            markers_taxa.loc[leaf.genome, "completeness"] = completenes
                            markers_taxa.loc[leaf.genome, "reference"] = genus


markers_taxa.to_csv(snakemake.output[0], index=True, sep="\t")
