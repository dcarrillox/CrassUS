import pandas as pd
import numpy as np
import os, glob, multiprocessing

# envs: utils.yaml

# ----------------------------------------------------------------------------
# read "aggregated_signals.txt" and init a dict to the store different signals
aggregated_df = pd.read_csv(snakemake.input[0], sep="\t", header=0, index_col=0)
#file = "/home/danielc/projects/CrassUS/results/ancient/aggregated_signals.txt"
#aggregated_df = pd.read_csv(file, sep="\t", header=0, index_col=0)

contigs_signals = {contig: {
                            "phylogenies": {
                                        "family":list(),
                                        "subfamily":list(),
                                        "genus":list()
                            },
                            "shared_prots": {
                                        "family":list(),
                                        "genus":list()
                            },
                            "ani": {
                                        "genus":list(),
                                        "species":list()
                            }
                        }
                        for contig in aggregated_df.index
                    }



# ---------------------------------
# Parse markers pylogeny annotation
def get_markers_annot(genome, markers, df):
    exclude = ["not_found", "too_short"]
    families    = list(set([df.loc[genome, f"family_{marker}"] for marker in markers if df.loc[genome, f"family_{marker}"] not in exclude]))
    subfamilies = list(set([df.loc[genome, f"subfamily_{marker}"] for marker in markers if df.loc[genome, f"subfamily_{marker}"] not in exclude]))
    genera      = list(set([df.loc[genome, f"genus_{marker}"] for marker in markers if df.loc[genome, f"genus_{marker}"] not in exclude]))

    family = list()
    subfam = list()
    genus  = list()

    if families:
        # check family by markers
        if len(families) == 1:
            family = families if families[0] != "unknown" else ["unknown"]
        else:
            if "unknown" in families:
                families.remove("unknown")
            family = families

        # check subfamily by markers
        if len(subfamilies) == 1:
            subfam = subfamilies if subfamilies[0] != "unknown" else ["unknown"]
        else:
            if "unknown" in subfamilies:
                subfamilies.remove("unknown")
            subfam = subfamilies

        # check genus by markers
        if len(genera) == 1:
            genus = genera if genera[0] != "unknown" else ["unknown"]
        else:
            if "unknown" in genera:
                genera.remove("unknown")
            genus = genera


    return family, subfam, genus

markers = sorted([marker for marker in snakemake.config["phylogenies"]["markers"] if snakemake.config["phylogenies"]["markers"][marker]], reverse=True)
#markers = ["TerL", "portal", "MCP"]
for genome in aggregated_df.index:
    markers_annot = get_markers_annot(genome, markers, aggregated_df)
    contigs_signals[genome]["phylogenies"]["family"] = markers_annot[0]
    contigs_signals[genome]["phylogenies"]["subfamily"] = markers_annot[1]
    contigs_signals[genome]["phylogenies"]["genus"] = markers_annot[2]



# ----------------------------
# Parse protein shared content
def get_protshared_annot(genome, df):
    family_cutoff = float(snakemake.config["shared_prots_cutoffs"]["family"])
    genus_cutoff = float(snakemake.config["shared_prots_cutoffs"]["genus"])

    families = list()
    genera = list()

    value_shared = df.loc[genome, "shared_prot_ref"]
    if not pd.isnull(value_shared):
        if round(value_shared) >= family_cutoff:
            families = df.loc[genome, "prot_most_similar_ref_family"].split(",")

            if value_shared >= genus_cutoff:
                genera = df.loc[genome, "prot_most_similar_ref_genus"].split(",")

    return families, genera

for genome in aggregated_df.index:
    protshared_annot = get_protshared_annot(genome, aggregated_df)
    contigs_signals[genome]["shared_prots"]["family"] = protshared_annot[0]
    contigs_signals[genome]["shared_prots"]["genus"] = protshared_annot[1]



# ---------------------
# Parse ani assignments
def get_ani_annot(genome, df):

    genus = str()
    species = str()

    # Genus
    ani_qcov_value = df.loc[genome, "ani_qcov"]
    ani_pid_value  = df.loc[genome, "ani_pid"]
    if not pd.isnull(ani_qcov_value) and not pd.isnull(ani_pid_value):
        if round(ani_qcov_value) >= 50:
            genus = df.loc[genome, "genus_most_similar_ref_genome"]
        else:
            genus = df.loc[genome, "ani_genus"]

        # Species
        if round(ani_qcov_value) >= 85 and round(ani_pid_value) >= 95:
            species = df.loc[genome, "species_most_similar_ref_genome"]
        else:
            species = df.loc[genome, "ani_species"]

    return genus, species

for genome in aggregated_df.index:
    ani_annot = get_ani_annot(genome, aggregated_df)
    contigs_signals[genome]["ani"]["genus"] = ani_annot[0]
    contigs_signals[genome]["ani"]["species"] = ani_annot[1]



# ------------------------
# init the final DataFrame
columns = ["crassus_id", "contig", "sample", "length", "len/taxa_len", "ref_taxa", "DTR",
           "family", "subfamily", "genus", "species",
           "evidence_family", "evidence_genus", "notes"]
final_df = pd.DataFrame(columns=columns)
final_df["crassus_id"] = aggregated_df.index
final_df.set_index("crassus_id", inplace=True)


# read tables with the contigs identifiers relation
def parse_id_table(ids_table):
    ids_table_df = pd.read_csv(ids_table, header=0, sep="\t", low_memory=False)
    return ids_table_df

pool = multiprocessing.Pool(processes=5)
ids_tables = glob.glob(f"{snakemake.params.ids_dir}/*.table")
ids_tables_dfs = pool.map(parse_id_table, ids_tables)
pool.close()
pool.join()

ids_df = pd.concat(ids_tables_dfs)
# remove rows with NaN, ie. contigs without new_id because they are short and were not analyzed
ids_df.dropna(inplace=True)
ids_df.set_index("new_id", inplace=True)


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

# add to final_df
for genome in final_df.index:
    # get raw id
    final_df.loc[genome, "contig"] = ids_df.loc[genome, "old_id"]
    # get sample
    final_df.loc[genome, "sample"] = genome.split("_")[0]
    # get length
    final_df.loc[genome, "length"] = int(genome.split("_")[-1])
    # get DTR
    final_df.loc[genome, "DTR"] = True if genome in dtr_genomes else False


# read taxas_average_length file, store to dict
lines = [line.strip().split("\t") for line in open(snakemake.params.lengths).readlines()]
taxas_lengths = {line[0]:float(line[1]) for line in lines}


# ----------------------------------------------------------------------
# iterate the assignments with the three signals and make the final call
def assess_family(contigs_signals, genome):
    # Family
    family = str()
    family_evidence = str()
    note = str()

    marker = contigs_signals[genome]["phylogenies"]["family"]
    protshared = contigs_signals[genome]["shared_prots"]["family"]

    # one or none family was predicted by markers
    if len(marker) <= 1:
        # markers and prot_shared are identical
        if marker == protshared:
            if marker and protshared: # make sure there is an annotation
                family = marker[0]
                family_evidence = "phylogenies, shared_prots"
            else:
                family = "NA"
        # markers and prot_shared are not identical
        else:
            # unknown by marker, known by protshared (only one family predicted)
            if (marker == ["unknown"] or not marker) and len(protshared) == 1:
                family = protshared[0]
                family_evidence = "shared_prots"
            # kwnon by marker, unknown by protshared
            elif (marker != ["unknown"] and marker) and not protshared:
                family = marker[0]
                family_evidence = "phylogenies"

    else:
        note = "multiple families in phylogenies"

    return family, family_evidence, note

def assess_subfamily(contigs_signals, genome):
    # Subfamily
    subfamily = str()
    note = str()

    marker = contigs_signals[genome]["phylogenies"]["subfamily"]

    # one or none family was predicted by markers
    if len(marker) <= 1:
        if marker:
            subfamily = marker[0]
    else:
        note = "multiple subfamilies in phylogenies"

    return subfamily, note

def assess_genus(contigs_signals, genome):
    # Genus
    genus = str()
    genus_evidence = str()
    note = str()

    marker = contigs_signals[genome]["phylogenies"]["genus"]
    protshared = contigs_signals[genome]["shared_prots"]["genus"]
    ani = contigs_signals[genome]["ani"]["genus"]

    # there is a ref genus assigned by ani
    genus_evidence = "ani"
    genus = ani
    if "genus__" not in ani:
        # check if it matches with protshared
        # if matches, assign that genus
        if ani in protshared:
            genus_evidence = "ani, shared_prots"
            if ani in marker:
                genus_evidence = "ani, shared_prots, phylogenies"
    else:
        if not protshared:
            genus_evidence = "ani, shared_prots"
            if not marker:
                genus_evidence = "ani, shared_prots, phylogenies"


    return genus, genus_evidence

def assess_species(contigs_signals, genome):
    # Species
    species = str()
    note = str()

    ani = contigs_signals[genome]["ani"]["species"]

    species = ani

    return species, note

def assess_completeness(df, taxas_lengths, genome):

    species = df.loc[genome, "species"]
    genus = df.loc[genome, "genus"]
    subfamily = df.loc[genome, "subfamily"]
    family = df.loc[genome, "family"]

    # print(genome)
    # print(species)
    # print(family)

    if not pd.isnull(species) and "species__" not in species:
        deepest = species
    elif not pd.isnull(genus) and "genus__" not in genus:
        deepest = genus
    elif not pd.isnull(subfamily) and subfamily != "unknown":
        deepest = subfamily
    elif not pd.isnull(family) and family != "unknown":
        deepest = family
    else:
        deepest = str()

    # calculate the completeness proxy
    if deepest:
        proxy = round(df.loc[genome, "length"]/taxas_lengths[deepest], 3)

        df.loc[genome, "len/taxa_len"] = proxy
        df.loc[genome, "ref_taxa"] = deepest



for genome in final_df.index:

    notes = list()

    # Family
    family, family_evidence, note = assess_family(contigs_signals, genome)
    final_df.loc[genome, "family"] = family
    final_df.loc[genome, "evidence_family"] = family_evidence
    if note:
        notes.append(note)

    # Subfamily
    subfamily, note = assess_subfamily(contigs_signals, genome)
    final_df.loc[genome, "subfamily"] = subfamily
    if note:
        notes.append(note)

    # Genus
    genus, genus_evidence = assess_genus(contigs_signals, genome)
    final_df.loc[genome, "genus"] = genus
    final_df.loc[genome, "evidence_genus"] = genus_evidence

    # Species
    species, note = assess_species(contigs_signals, genome)
    final_df.loc[genome, "species"] = species
    if note:
        notes.append(note)


    # write notes if necessary
    if notes:
        final_df.loc[genome, "notes"] = "; ".join(notes)


    # write completeness
    assess_completeness(final_df, taxas_lengths, genome)



# Write to final output file
final_df.to_csv(snakemake.output[0], sep="\t", header=True, index=True)






















# for genome in contigs_signals:
#     print(genome, contigs_signals[genome])
#     print()






#
# def assess_genus(genome, df):
#     check = False
#     # check the prots shared are higher than the cutoff for the genus
#     if df.loc[genome, "shared_prot_ref"] >= float(snakemake.config["shared_prots_cutoffs"]["genus"]):
#         shared_genera = df.loc[genome, "prot_most_similar_ref_genus"].split(",")
#         # check the AF is higher than the cutoff
#         if df.loc[genome, "qcov"] >= 50:
#             # check if the genus called by AF agrees with one of the called by shared_prots
#             af_genus = df.loc[genome, "AF_most_similar_ref_genus"]
#             if af_genus in shared_genera:
#                 genus = af_genus
#                 check = True
#     if not check:
#         genus = df.loc[genome, "AF_genus"]
#
#     return genus
#
#
#

#
# dtr_genomes = set(dtr_genomes)
#
# # read average lengths for each rank, store to dict
# lines = [line.strip().split("\t") for line in open(snakemake.params.lengths).readlines()]
# taxas_lengths = {line[0]:float(line[1]) for line in lines}
#
# ##
# # read aggregated results
# aggregated_df = pd.read_csv(snakemake.input[0], sep="\t", header=0, index_col=0)
#
# # get markers from the table itself
# markers = [column.split("family_")[1] for column in aggregated_df.columns if column.startswith("family_")]
#
# # init the dict to store the information when iterating the aggregated table: k=genome  v=columns
# to_write = {genome :{"genome": genome,
#                      "len/taxa_len":"",
#                      "ref_taxa":str(),
#                      "DTR":"no",
#                      "family":str(),
#                      "subfamily":str(),
#                      "genus":str(),
#                      "species":str(),
#                      "family_evidence":list(),
#                      "notes":list()
#                     }
#             for genome in aggregated_df.index}
#
#
#
#
#
# # iterate the genomes while assessing their final taxonomy
# for genome in aggregated_df.index:
#     genome_length = float(genome.split("_")[-1])
#     # get family and subfamily from markers
#     marker_fam, marker_subfam, marker_genus = get_markers_annot(genome, markers, aggregated_df)
#
#     # check if shared prot validates the marker assignment
#     # to_write[genome]["family"] = marker_fam
#     # to_write[genome]["subfamily"] = marker_subfam
#
#     # check family given by shared_prots
#     family_shared = family_by_shared(genome, aggregated_df)
#     # compare the family/ies to the marker one. If the later is unknown, then we can call the shared one
#     if family_shared:
#         if marker_fam not in ["", "unknown"]:
#             if marker_fam not in family_shared:
#                 to_write[genome]["notes"].append("shared prots with family below 20%")
#             else:
#                 to_write[genome]["family"] = marker_fam
#                 to_write[genome]["subfamily"] = marker_subfam
#                 to_write[genome]["family_evidence"].append("phylogenies")
#                 to_write[genome]["family_evidence"].append("shared prots")
#         else:
#             if len(family_shared) == 1:
#                 to_write[genome]["family"] = family_shared[0]
#                 to_write[genome]["subfamily"] = marker_subfam
#                 to_write[genome]["family_evidence"].append("shared prots")
#     else:
#         to_write[genome]["family"] = marker_fam
#         to_write[genome]["subfamily"] = marker_subfam
#         to_write[genome]["family_evidence"].append("phylogenies")
#         if marker_fam not in ["", "unknown"]:
#             to_write[genome]["notes"].append("shared prots with family below 20%")
#
#
#
#
#     # get genus & species
#     genus = assess_genus(genome, aggregated_df)
#     to_write[genome]["genus"] = genus
#     to_write[genome]["species"] = aggregated_df.loc[genome, "ani_species"]
#
#     # check if DTR
#     if genome in dtr_genomes:
#         to_write[genome]["DTR"] = "yes"
#
#     # assess completeness based on the lowest rank assigned
#     compl = ""
#     if not genus.startswith("genus__"):
#         to_write[genome]["ref_taxa"] = genus
#         compl = round(genome_length/taxas_lengths[genus] , 2)
#     else:
#         if to_write[genome]["subfamily"] not in  ["", "unknown"]:
#             to_write[genome]["ref_taxa"] = to_write[genome]["subfamily"]
#             compl = round(genome_length/taxas_lengths[to_write[genome]["subfamily"]] , 5)
#         elif to_write[genome]["family"] not in  ["", "unknown"]:
#             to_write[genome]["ref_taxa"] = to_write[genome]["family"]
#             compl = round(genome_length/taxas_lengths[to_write[genome]["family"]] , 5)
#     to_write[genome]["len/taxa_len"] = compl
#
#     # join notes & family_evidence
#     to_write[genome]["notes"] = ";".join(to_write[genome]["notes"])
#     to_write[genome]["family_evidence"] = ";".join(to_write[genome]["family_evidence"])
#
#
# to_write_df = pd.DataFrame.from_dict(to_write, orient="index")
# to_write_df.to_csv(snakemake.output[0], index=False, sep="\t")
