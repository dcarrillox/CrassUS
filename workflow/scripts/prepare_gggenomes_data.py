import pandas as pd
import os, glob


os.makedirs(snakemake.output[0], exist_ok=True)


# check if reference
def is_ref(tname):
    return tname in crass_reference

# Group by query genome and select targets
def filter_func(query_df):
    sorted_df = query_df.sort_values(by='qcov', ascending=False)
    qname = sorted_df.iloc[0]['qname']
    tnames = sorted_df['tname'].to_list()

    # Check if first one is reference
    first_one = tnames[0]
    if is_ref(first_one):
        df = pd.DataFrame.from_dict( {'qname': [qname],  'tname':[first_one], 'group': ['ref']} )
        return df
    else:
        try:
            first_ref = tnames[[is_ref(tname) for tname in tnames].index(True)]
            df = pd.DataFrame.from_dict({'qname': [qname, qname],  'tname':[first_one, first_ref], 'group': ['noref','ref']})
            return df
        except ValueError:
            return pd.DataFrame.from_dict({'qname': [qname],  'tname':[''], 'group': ['']}) # no ANI matches

def make_straight_alignments(key_df):
    if key_df["tstart_tmp"].iloc[0] > key_df["tend_tmp"].iloc[0]:
        key_df.rename({'tstart_tmp': 'tend', 'tend_tmp': 'tstart'}, axis=1, inplace=True)
    else:
        key_df.rename({'tstart_tmp': 'tstart', 'tend_tmp': 'tend'}, axis=1, inplace=True)
    return key_df

def process_genome_table(genome_id, crass_reference, func_annot, target_df):
    if genome in crass_reference:
        file = f"{snakemake.params.ref_genome_tables_dir}/{genome}.table"
    else:
        file = glob.glob(f"{snakemake.params.genome_tables_dir}/{genome}_tbl*.table")[0]


    df = pd.read_csv(file, sep="\t", header=0, index_col=0)
    for protein in df.index:
        if df.loc[protein, "strand"] == "-":
            start, end = df.loc[protein, "start"], df.loc[protein, "end"]
            df.loc[protein, "start"] = end
            df.loc[protein, "end"]  = start

        if pd.isnull(df.loc[protein, "yutin"]):
            df.loc[protein, "yutin"] = "NA"
        else:
            if df.loc[protein, "yutin"] not in func_annot:
                df.loc[protein, "yutin"] = "Others"

    return df




# ------------------------------------------------------------------------------------
# read and filter anicalc results to keep only non-reference queries and non-self hits
anicalc_df = pd.read_csv(snakemake.input.ani[0], sep="\t", header=0)

# Load crass referemces
crass_reference = {line.split("\t")[0] for line in open(snakemake.params.taxonomy).readlines()}

# Remove reference queries
anicalc_df = anicalc_df[~anicalc_df.qname.isin(crass_reference)]



# ------------------------
# Group anicalc results by query, non-ref contig. Keep first hit. If this is a
# reference genome, that's it. Otherwise, keep the first genome + the closest
# reference genome.
res = anicalc_df.groupby('qname').apply(filter_func)
res['key'] = res['qname'] + '@' + res['tname']
res = res[['key', 'group']]



# -----------------------------------------------
# Read blast file. Add column "key" so I can filter with it
names = ["qname", "tname", "perc_sim", "aln_len", "mismatch", "gap", "qstart", "qend", "tstart_tmp", "tend_tmp", "evalue", "std", "qlen", "slen"]
blast_df = pd.read_csv(snakemake.input.blast, sep="\t", names=names)
# remove alignments shorter than the cutoff
min_aln_len = int(snakemake.config["plot"]["aln_min_len"])
blast_df = blast_df[blast_df["aln_len"] >= min_aln_len]
blast_df['key'] = blast_df['qname'] + '@' + blast_df['tname']


# filter df based on the "keys" obtained from anicalc results
ani_key_set = set(res.key.values)
filtered_df = blast_df[blast_df.key.isin(ani_key_set)]

# join back the groups
filtered_df = filtered_df.merge(res, on = 'key', how = 'left')
query_groups = filtered_df.groupby('qname')

func_annot = ["TerL", "MCP", "portal", "primase", "DNApB", "PolA", "Tstab", "Ttub", "tail_fib"]



# Write to output
for query, query_df in query_groups:
    query_target_pair = query_groups.get_group(query)

    # group by key and make alignments straight for each of them
    key_groups = query_df.groupby("key").apply(make_straight_alignments)
    key_groups.to_csv(f"{snakemake.output[0]}/{query}.blast", sep="\t", index=False)

    # get genome tables for the genomes involved
    genomes = [query] + list(set(query_target_pair["tname"]))
    tables_df = pd.DataFrame()
    for genome in genomes:
        df = process_genome_table(genome, crass_reference, func_annot, tables_df)
        tables_df = tables_df.append(df, ignore_index=True)
    tables_df.to_csv(f"{snakemake.output[0]}/{query}.annot", sep="\t")
