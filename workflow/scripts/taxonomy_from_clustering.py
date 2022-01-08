import pandas as pd

# env: utils.yaml

# read shared content matrix. It contains values for crassus genomes only.
matrix_shared = pd.read_csv(snakemake.input.matrix_shared, sep="\t", header=0, index_col=0)


# read reference taxonomy
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()[1:]]
for line in lines:
    crass_taxonomy[line[0]] = {"family":line[1], "subfamily":line[2], "genus":line[3]}

to_write = list()
# find out which (complete) genomes share >70% with any other genome
for qgenome in matrix_shared.columns[4:]:
    if qgenome not in crass_taxonomy:
        to_add = [qgenome]
        # sort the df by the values for this query genome
        matrix_shared = matrix_shared.sort_values(qgenome, ascending=False)

        # get the most similar ref genome
        ref_family_list = list()
        ref_subfam_list = list()
        ref_genus_list  = list()
        max_shared = 0
        for tgenome in matrix_shared.index:
            if tgenome in crass_taxonomy:
                value = round(matrix_shared.loc[tgenome,qgenome] * 100,2)
                if value >= max_shared:
                    max_shared = value
                    ref_family_list.append(crass_taxonomy[tgenome]["family"])
                    ref_subfam_list.append(crass_taxonomy[tgenome]["subfamily"])
                    ref_genus_list.append(crass_taxonomy[tgenome]["genus"])

        if ref_family_list and ref_genus_list:
            if max_shared != 0:
                to_add.append(max_shared)
                ref_family = ",".join(sorted(list(set(ref_family_list))))
                to_add.append(ref_family)
                ref_subfam = ",".join(sorted(list(set(ref_subfam_list))))
                to_add.append(ref_subfam)
                ref_genus  = ",".join(sorted(list(set(ref_genus_list))))
                to_add.append(ref_genus)
            else:
                to_add.append("")
                to_add.append("")
                to_add.append("")
        else:
            to_add.append("")
            to_add.append("")
            to_add.append("")

        to_write.append(to_add)


# create df and write to file
columns = ["genome", "ref_shared_prot", "most_similar_family", "most_similar_subfamily", "most_similar_genus"]
to_write_df = pd.DataFrame(to_write, columns=columns)
to_write_df.to_csv(snakemake.output[0], index=False, sep="\t")
