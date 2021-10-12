import pandas as pd

# env: utils.yaml

# read shared content matrix. It contains values for crassus genomes only.
matrix_shared = pd.read_csv(snakemake.input.matrix_shared, sep="\t", header=0, index_col=0)


# read reference taxonomy
reference_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()]
for line in lines:
    reference_taxonomy[line[0]] = {"family":line[2], "subfamily":line[3], "genus":line[4]}


# read completeness table to know which genomes were classified with markers at the genus level
markers_table = pd.read_csv(snakemake.input.markers_table[0], sep="\t", header=0, index_col=0)
markers = [marker for marker in snakemake.config["phylogenies"] if snakemake.config["phylogenies"][marker]]
# store in a dict the genomes that got genus marker classification. For this, iterate
# the rows and check the content of "genus_MARKER" columns
no_taxa_assigned = ["Not found",
                    "truncated",
                    "wrong strands, check func. annot",
                    "multiple copies",
                    "unknown"]
genus_marker = dict()
genus_marker_columns = [f"genus_{marker}" for marker in markers]
for genome in markers_table.index:
    genera = list()
    for column in genus_marker_columns:
        genera.append(markers_table.loc[genome, column])
    # remove not found and unknown. Also, remove _1 from monophyletic correction
    final_genus = list(set([genus.split("_")[0] for genus in genera if genus not in no_taxa_assigned]))
    if final_genus:
        # check that the different markers agree
        if len(final_genus) == 1:
            genus_marker[genome] = final_genus[0]
        else:
            print(f"{genome}, discrepancies at the genus level -> {final_genus}")


to_write = list()
# get complete AND classified genomes according to markers (completeness) table
complete_genomes = markers_table[(markers_table.completeness >= 90) & (markers_table.highest_taxa.notnull())].index.tolist()

# shared cuttof for genus delimitation
shared_genus_cutoff = float(snakemake.config["shared_genus_cutoff"])

# find out which (complete) genomes share >70% with any other genome
for query_genome in matrix_shared:
    # query_genome is complete
    if query_genome in complete_genomes:
        # sort the df by the values for this query genome
        matrix_shared = matrix_shared.sort_values(query_genome, ascending=False)
        # get genomes >70%
        shared07_genomes = matrix_shared[matrix_shared[query_genome] >= shared_genus_cutoff].index.tolist()
        # there are similar genomes
        if shared07_genomes:
            shared_genera = list()
            for shared07_genome in shared07_genomes:
                # check if the target genome is a reference crassphage
                if shared07_genome in reference_taxonomy:
                    shared_genera.append(reference_taxonomy[shared07_genome]["genus"])
                # otherwise, check if the genome obtained genus marker
                elif shared07_genome in genus_marker:
                    shared_genera.append(genus_marker[shared07_genome])


            genera = list(set(shared_genera))
            n_genera = len(genera)

            # if not similar genus was found, set it as "new_genus"
            if not genera:
                genera = ["new_genus"]

            # get the most similar genome in the reference, so it has genus
            most_similar_genus = "new_genus"
            for shared07_genome in shared07_genomes:
                if shared07_genome in reference_taxonomy:
                    most_similar_genus = reference_taxonomy[shared07_genome]["genus"]
                    break

            to_write.append([query_genome, "yes", len(shared07_genomes), n_genera, ",".join(genera), most_similar_genus, ",".join(shared07_genomes)])



        # there are NOT similar genomes
        else:
            to_write.append([query_genome, "yes", 0, 0, "new_genus", "new_genus",  ",".join(shared07_genomes)])

    # query_genome is not complete
    else:
        to_write.append([query_genome, "no", "", "", "", "", ""])




# create df and write to file
columns = ["genome", "complete", "n_shared70", "n_genera", "genera", "most_similar_genus", "most_similar_genomes"]
to_write_df = pd.DataFrame(to_write, columns=columns)
#print(to_write_df)
to_write_df.to_csv(snakemake.output[0], index=False, sep="\t")
