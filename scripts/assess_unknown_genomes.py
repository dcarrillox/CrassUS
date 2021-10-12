from ete3 import Tree
import pandas as pd
import os, glob

# env: phylogenies.yaml

# read taxa assignments and get the unknown genomes
markers_taxa = pd.read_csv(snakemake.input.taxa_assessments[0], sep="\t", header=0, index_col=0)
unknown = markers_taxa[markers_taxa.highest_taxa == "unknown"]

# for each of the unknown genomes, read its genome table and store the number of proteins
geno_tables_dir = snakemake.params.geno_tables_dir
genomes_nprots = dict()
for genome in unknown.index:
    geno_table_file = glob.glob(f"{geno_tables_dir}/{genome}_tbl-*.table")[0]
    nprots = len(open(geno_table_file).readlines())-1
    genomes_nprots[genome] = nprots



# read crass taxonomy, store only the family
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
for line in lines:
    crass_taxonomy[line[0]] = line[2]

taxa_to_family = dict()
for line in lines:
    taxa_to_family[line[2]] = line[2]
    taxa_to_family[line[3]] = line[2]
    taxa_to_family[line[4]] = line[2]

# get crassus taxonomy, store only the family
crassus_taxonomy = dict()
crassus_classified = markers_taxa[(markers_taxa.highest_taxa != "unknown") & (markers_taxa.highest_taxa.notnull())]
for genome in crassus_classified.index:
    family = taxa_to_family[crassus_classified.loc[genome, "highest_taxa"]]
crassus_taxonomy = {genome:family for genome in crassus_classified.index}


# read sharing percentages
sharing_matrix = pd.read_csv(snakemake.input.sharing_percentages, sep="\t", header=0, index_col=0)
# iterate the unknown genomes and check with which genome they share the highest amount of proteins
to_write_df = pd.DataFrame(columns=["nprots", "most_sim_genome", "most_sim_fam", "sharing",
                                    "most_sim_classfd_genome","most_sim_classfd_fam", "sharing_classfd"],
                           index=unknown.index )
for genome in unknown.index:
    to_write_df.loc[genome, "nprots"] = genomes_nprots[genome]
    # sort sharing matrix by the unknown genome
    sharing_matrix = sharing_matrix.sort_values(genome, ascending=False)
    # most similar genome
    most_sim_genome = sharing_matrix[genome].index[0]
    to_write_df.loc[genome, "most_sim_genome"] = most_sim_genome
    to_write_df.loc[genome, "sharing"] = sharing_matrix[genome][0]
    if most_sim_genome in crass_taxonomy:
        to_write_df.loc[genome, "most_sim_fam"] = crass_taxonomy[most_sim_genome]
    elif most_sim_genome in crassus_taxonomy:
        to_write_df.loc[genome, "most_sim_fam"] = crassus_taxonomy[most_sim_genome]

    # get the most similar genome classified
    for target_genome in sharing_matrix.index:
        if target_genome in crass_taxonomy:
            to_write_df.loc[genome, "most_sim_classfd_genome"] = target_genome
            to_write_df.loc[genome, "most_sim_classfd_fam"] = crass_taxonomy[target_genome]
            to_write_df.loc[genome, "sharing_classfd"] = sharing_matrix.loc[target_genome, genome]
            break

        elif most_sim_genome in crassus_taxonomy:
            to_write_df.loc[genome, "most_sim_classfd_genome"] = target_genome
            to_write_df.loc[genome, "most_sim_classfd_fam"] = crassus_taxonomy[target_genome]
            to_write_df.loc[genome, "sharing_classfd"] = sharing_matrix.loc[target_genome, genome]
            break


# prepare the df to write it, but don't do it yet. First remove unknown contigs from the trees
to_write_df = to_write_df.sort_values("nprots", ascending=False)


# read trees and remove unclassified genomes
for marker_tree in snakemake.input.markers_trees:
    # read tree
    t = Tree(marker_tree, format=1)
    # assign taxonomy for the reference crAssphages
    for leaf in t.iter_leaves():
        # check if the leaf comes from the reference set
        genome = leaf.name.split("|")[0]
        if genome in crass_taxonomy:
            leaf.add_features(family=crass_taxonomy[genome])
        else:
            if genome in to_write_df.index:
                leaf.add_features(family="unknown")
            else:
                leaf.add_features(family=crassus_taxonomy[genome])

    # find the LCA of the two outgroup species
    outgs_leaves = t.search_nodes(family="outgroup")
    outgs_lca = t.get_common_ancestor(outgs_leaves)
    # reroot the tree
    t.set_outgroup(outgs_lca)

    lcas = t.get_monophyletic(values=["unknown"], target_attr="family")
    if lcas:
        for lca in lcas:
            lca.detach()

    t.write(format=1, outfile=marker_tree.replace(".nwk", "_classified.nwk"))


to_write_df.to_csv(snakemake.output[0], index=True, sep="\t")
