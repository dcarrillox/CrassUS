import pandas as pd
import os

# env: utils.yaml

# start by parsing the DTR results
# input is the mock .dtr_blast_done file, the actual file is .dtr_blast
# only DTR genomes will be in the dictionary
genomes_completeness = dict()
for dtr_file in snakemake.input.dtr_blast_done:
    # grab the actual blast file
    blast_file = dtr_file.replace("_done", "")
    # get the genome, it will be the key of the dictionary
    genome = os.path.basename(blast_file).replace(".dtr_blast", "")
    # check it is >20Kb
    genome_length = int(genome.split("_")[-1])
    if genome_length >= 20000:
        lines = [line.strip().split("\t") for line in open(blast_file).readlines()]
        # discard self-hits
        lines = [line for line in lines if line[0] != line[1]]
        # check if query hit
        for line in lines:
            if "|start" in line[0] and "|end" in line[1]:
                similarity = float(line[2])
                aln_length = int(line[3])
                query_start = int(line[6])
                if similarity >= 98 and aln_length >= 20 and query_start <= 50:
                    genomes_completeness[genome] = [100, "DTR"]



# for the non DTR, read average lengths for each rank, store to dict
lines = [line.strip().split("\t") for line in open(snakemake.params.lengths).readlines()]
taxas_lengths = {line[0]:float(line[1]) for line in lines}

# read taxa assignment from the trees
df = pd.read_csv(snakemake.input.taxa_markers[0], sep="\t", index_col = 0)
# get the markers from config
markers = [marker for marker in snakemake.config["phylogenies"] if snakemake.config["phylogenies"][marker]]
# get the genomes
genomes = df.index.to_list()

ranks = ["genus", "subfamily", "family"]

for genome in genomes:
    for rank in ranks:
        # genera & subfamily were adjusted for some ref genomes to make them monophyletic,
        # split by "_" to get the actual name, ie. Blohavirus instead of Blohavirus_1
        assigned_taxas = [df.loc[genome, f"{rank}_{marker}"].split("_")[0] for marker in markers]
        # remove redundancy
        assigned_taxas = list(set(assigned_taxas))
        # remove "not found" and "unknown"
        if "Not found" in assigned_taxas:
            assigned_taxas.remove("Not found")
        if "unknown" in assigned_taxas:
            assigned_taxas.remove("unknown")
        # if at the end there is a taxa assigned, use it
        # check if it is more than one tho
        if assigned_taxas:

            if len(assigned_taxas) == 1:
                # this is the assigned taxonomy by the marker tree
                # get the average length for it
                ref_len = taxas_lengths[assigned_taxas[0]]
                genome_len = int(genome.split("_")[-1])
                if genome not in genomes_completeness:
                    compl_value = round((genome_len*100)/ref_len, 2)
                    #genomes_completeness[genome] = float(genome_len/ref_len)
                    genomes_completeness[genome] = [compl_value, assigned_taxas[0]]
                break
            else:
                print(genome, f"discrepancies between the different markers at the {rank} level.")


#
df.insert(0, 'reference',"")
df.insert(0, 'completeness',"")

for genome in genomes:
    if genome not in genomes_completeness:
        df.loc[genome, "completeness"] = "NaN"
        df.loc[genome, "reference"] = "NaN"
    else:
        df.loc[genome, "completeness"] = genomes_completeness[genome][0]
        df.loc[genome, "reference"] = genomes_completeness[genome][1]


print(df)
df.to_csv(snakemake.output[0], sep="\t", index=True)


#with open(snakemake.output[0], "w") as fout:
#os.system(f"touch {snakemake.output[0]}")
