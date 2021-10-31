import pandas as pd
import os

# envs: utils.yaml

# read DTR blast, store results to list
dtr_genomes = list()
for dtr_file in snakemake.input.dtr_blast_done:
    # grab the actual blast file
    blast_file = dtr_file.replace("_done", "")
    # get the genome, it will be the key of the dictionary
    genome = os.path.basename(blast_file).replace(".dtr_blast", "")
    # check it is >20Kb
    genome_length = int(genome.split("_")[-1])
    if genome_length >= 10000:
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
                    dtr_genomes.append(genome)

dtr_genomes = set(dtr_genomes)

# read average lengths for each rank, store to dict
lines = [line.strip().split("\t") for line in open(snakemake.params.lengths).readlines()]
taxas_lengths = {line[0]:float(line[1]) for line in lines}

##
# read aggregated results
aggregated_df = pd.read_csv(snakemake.input[0], sep="\t", header=0, index_col=0)

# get markers from the table itself
markers = [column.split("family_")[1] for column in aggregated_df.columns if column.startswith("family_")]

# init the dict to store the information when iterating the aggregated table: k=genome  v=columns
to_write = {genome :{"genome": genome,
                     "len/taxa_len":"",
                     "ref_taxa":str(),
                     "DTR":"no",
                     "family":str(),
                     "subfamily":str(),
                     "genus":str(),
                     "species":str(),
                     "notes":list()
                    }
            for genome in aggregated_df.index}


def get_markers_annot(genome, markers, df):
    exclude = ["not_found", "too_short", "unknown"]
    family    = list(set([df.loc[genome, f"family_{marker}"] for marker in markers if df.loc[genome, f"family_{marker}"] not in exclude]))
    subfamily = list(set([df.loc[genome, f"subfamily_{marker}"] for marker in markers if df.loc[genome, f"subfamily_{marker}"] not in exclude]))
    if len(family) == 1 and len(subfamily) == 1:
        return family[0], subfamily[0]
    else:
        return "", ""

def assess_genus(genome, df):
    check = False
    # check the prots shared are higher than the cutoff for the genus
    if df.loc[genome, "prot_ref"] >= float(snakemake.config["shared_prots_cutoffs"]["genus"]):
        shared_genera = df.loc[genome, "prot_most_similar_ref_genus"].split(",")
        # check the AF is higher than the cutoff
        if df.loc[genome, "AF"] >= 50:
            # check if the genus called by AF agrees with one of the called by shared_prots
            af_genus = df.loc[genome, "AF_most_similar_ref_genus"]
            if af_genus in shared_genera:
                genus = af_genus
                check = True
    if not check:
        genus = df.loc[genome, "AF_genus"]

    return genus



# iterate the genomes while assessing their final taxonomy
for genome in aggregated_df.index:
    genome_length = float(genome.split("_")[-1])
    # get family and subfamily from markers
    marker_fam, marker_subfam = get_markers_annot(genome, markers, aggregated_df)

    # check if shared prot validates the marker assignment
    to_write[genome]["family"] = marker_fam
    to_write[genome]["subfamily"] = marker_subfam
    if marker_fam != "":
        if aggregated_df.loc[genome, "prot_ref"] >= 20 and marker_fam in aggregated_df.loc[genome, "prot_most_similar_ref_family"].split(","):
            pass
        else:
            to_write[genome]["notes"].append("shared prots with family below 20%")

    # get genus & species
    genus = assess_genus(genome, aggregated_df)
    to_write[genome]["genus"] = genus
    to_write[genome]["species"] = aggregated_df.loc[genome, "ani_species"]

    # check if DTR
    if genome in dtr_genomes:
        to_write[genome]["DTR"] = "yes"

    # assess completeness based on the lowest rank assigned
    compl = ""
    if not genus.startswith("genus__"):
        to_write[genome]["ref_taxa"] = genus
        compl = round(genome_length/taxas_lengths[genus] , 2)
    else:
        if marker_subfam != "":
            to_write[genome]["ref_taxa"] = marker_subfam
            compl = round(genome_length/taxas_lengths[marker_subfam] , 2)
        elif marker_fam != "":
            to_write[genome]["ref_taxa"] = marker_fam
            compl = round(genome_length/taxas_lengths[marker_fam] , 2)
    to_write[genome]["len/taxa_len"] = compl

    # join notes
    to_write[genome]["notes"] = ";".join(to_write[genome]["notes"])


to_write_df = pd.DataFrame.from_dict(to_write, orient="index")
to_write_df.to_csv("/home/danielc/projects/crAssUS/test.txt", index=False, sep="\t")
