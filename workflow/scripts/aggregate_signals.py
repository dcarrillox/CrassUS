import pandas as pd

# envs: utils.yaml


# read markers table, it will be the template where to add the rest of signals
df = pd.read_csv(snakemake.input.phylogenies[0], header=0, sep="\t", index_col=0)


# read shared_content
shared_df = pd.read_csv(snakemake.input.shared_prot, header=0, sep="\t", index_col=0)
for genome in shared_df.index:
    df.loc[genome, "prot_ref"] = shared_df.loc[genome, "ref_shared_prot"]
    df.loc[genome, "prot_most_similar_ref_family"] = shared_df.loc[genome, "most_similar_family"]
    df.loc[genome, "prot_most_similar_ref_genus"] = shared_df.loc[genome, "most_similar_genus"]


# read ani_clusters assignments
aniclust_df = pd.read_csv(snakemake.input.ani_cluster, header=0, sep="\t", index_col=0)
for genome in aniclust_df.index:
    df.loc[genome, "AF"] = aniclust_df.loc[genome, "AF"]
    df.loc[genome, "AF_genus"] = aniclust_df.loc[genome, "genus"]
    df.loc[genome, "AF_most_similar_ref_genus"] = aniclust_df.loc[genome, "most_similar_ref"]
    df.loc[genome, "ani_species"] = aniclust_df.loc[genome, "species"]

df.to_csv(snakemake.output[0], sep="\t", index=True)
