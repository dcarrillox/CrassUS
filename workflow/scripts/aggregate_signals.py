import pandas as pd

# envs: utils.yaml


# read markers table, it will be the template to add the rest of signals
df = pd.read_csv(snakemake.input.phylogenies[0], header=0, sep="\t", index_col=0)


# read coding prediction
coding_df = pd.read_csv(snakemake.input.coding, header=0, sep="\t")
coding_df = coding_df[coding_df["pick"] == True]
coding_df.set_index("contig", inplace=True)
for genome in coding_df.index:
    df.loc[genome, "coding"] = coding_df.loc[genome, "coding"]


# read shared_content
shared_df = pd.read_csv(snakemake.input.shared_prot, header=0, sep="\t", index_col=0)
for genome in shared_df.index:
    df.loc[genome, "shared_prot_ref"] = shared_df.loc[genome, "ref_shared_prot"]
    df.loc[genome, "prot_most_similar_ref_family"] = shared_df.loc[genome, "most_similar_family"]
    df.loc[genome, "prot_most_similar_ref_genus"] = shared_df.loc[genome, "most_similar_genus"]


# read ani_clusters assignments
aniclust_df = pd.read_csv(snakemake.input.ani_cluster, header=0, sep="\t", index_col=0)
for genome in aniclust_df.index:
    df.loc[genome, "ani_most_similar_genome"] = aniclust_df.loc[genome, "most_similar_genome"]
    df.loc[genome, "ani_most_similar_ref_genome"] = aniclust_df.loc[genome, "most_similar_ref_genome"]
    df.loc[genome, "genus_most_similar_ref_genome"] = aniclust_df.loc[genome, "ref_genus"]
    df.loc[genome, "species_most_similar_ref_genome"] = aniclust_df.loc[genome, "ref_species"]
    df.loc[genome, "ani_qcov"] = aniclust_df.loc[genome, "qcov"]
    pid = aniclust_df.loc[genome, "pid"]
    df.loc[genome, "ani_pid"] = pid if not pd.isnull(pid) else 0
    df.loc[genome, "ani_genus"] = aniclust_df.loc[genome, "genus_names"]
    df.loc[genome, "ani_species"] = aniclust_df.loc[genome, "species_names"]

df.to_csv(snakemake.output[0], sep="\t", index=True)
