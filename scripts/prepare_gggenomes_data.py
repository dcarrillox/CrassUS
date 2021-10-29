import pandas as pd
import os
import multiprocessing
from functools import partial
from collections import defaultdict
import time

os.makedirs(snakemake.output[0], exist_ok=True)

anicalc_df = pd.read_csv(snakemake.input.ani[0], sep="\t", header=0)

# Load crass referemces
crass_reference = {line.split("\t")[0] for line in open(snakemake.params.taxonomy).readlines()}

# Remove reference queryies
anicalc_df = anicalc_df[~anicalc_df.qname.isin(crass_reference)]



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


# df = pd.DataFrame.from_dict( {'qname': 'ee',  'tname':'test', 'group': 'ref'} )
res = anicalc_df.groupby('qname').apply(filter_func)
res['key'] = res['qname'] + '@' + res['tname']
res = res[['key', 'group']]

# load the blast file
names = ["qname", "tname", "perc_sim", "aln_len", "mismatch", "gap", "qstart", "qend", "tstart", "tend", "evalue", "std", "qlen", "slen"]
blast_df = pd.read_csv("/home/danielc/projects/crAssUS/results/7_ANI/0_species/all_genomes_ref_crassus_blast.tsv", sep="\t", names=names)
blast_df['key'] = blast_df['qname'] + '@' + blast_df['tname']

# again make key column
ani_key_set = set(res.key.values)

# subset keys
filtered_df = blast_df[blast_df.key.isin(ani_key_set)]

# join back the groups
filtered_df = filtered_df.merge(res, on = 'key', how = 'left')
query_groups = filtered_df.groupby('qname')

# Write to output
for x in query_groups.groups:
    query_target_pair = query_groups.get_group(x)
    query_target_pair.to_csv(f"{snakemake.output[0]}/{x}.txt", sep="\t", index=False)
    # print(query_target_pair)
    # print('----\n\n')
