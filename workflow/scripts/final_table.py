import pandas as pd
import numpy as np
import os, glob, multiprocessing

# envs: utils.yaml


# ----------------------------------------------------------------------------
# read "aggregated_signals.txt" and init a dict to the store different signals
aggregated_df = pd.read_csv(snakemake.input[0], sep="\t", header=0, index_col=0)

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

    # if genome == "Baboon19_894_46906":
    #     print(families)

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
        #     else:
        #         genera = ["unknown"]
        #
        # else:
        #     families = ["unknown"]
        #     genera = ["unknown"]

    # if genome == "Baboon19_140_109222":
    #     print(genome)
    #     print(families)
    #     print(genera)

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
           "evidence_family", "evidence_genus", "notes", "discard"]
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

    if genome == "Baboon19_708_52401":
        print(genome)
        print("marker:", marker)
        print("prot:", protshared)

    # one or none family were predicted by markers

    # print(genome)
    # print(marker)
    # print(protshared)
    # print()

    if len(marker) <= 1:
        # markers and prot_shared are identical
        if marker == protshared:
            # check there are assignments
            if marker and protshared:
                family = marker[0]
                family_evidence = "phylogenies, shared_prots"
            else:
                family = ""
                family_evidence = "phylogenies, shared_prots"

        # markers and prot_shared are not identical
        else:
            # not signal by markers
            if not marker:
                # check number of families predicted by shared prots
                if protshared:
                    if len(protshared) > 1:
                        family = "unknown"
                        family_evidence = "shared_prots"
                        note = "multiple families by shared prots"
                    if len(protshared) == 1:
                        family = protshared[0]
                        family_evidence = "shared_prots"
                    else:
                        family = ""

            else:
                # unknown by marker, known by protshared, one or more families
                if marker == ["unknown"]:
                    if len(protshared) > 1:
                        family = "unknown"
                        family_evidence = "phylogenies, shared_prots"
                        note = "multiple families by shared prots"
                    elif len(protshared) == 1:
                        family = protshared[0]
                        family_evidence = "shared_prots"
                    else:
                        family = "unknown"
                        family_evidence = "phylogenies, shared_prots"

                # kwnon by marker, unknown by protshared
                elif marker != ["unknown"] and not protshared:
                    family = marker[0]
                    family_evidence = "phylogenies"
                    note = "shared proteins below the family cutoff"
                # known by marker, multiple families by shared
                elif marker != ["unknown"] and len(protshared) > 1:
                    if marker[0] in protshared:
                        family = marker[0]
                        family_evidence = "phylogenies, shared_prots"
                        note = "multiple families by shared prots"


    else:
        note = "multiple families in phylogenies"

    # if genome == "Baboon19_1549_35692":
    #     print(family)

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

    # if genome == "Baboon19_140_109222":
    #     print(genome)
    #     print(marker)
    #     print(protshared)
    #     print(ani)

    # print(genome)
    # print(marker)
    # print(protshared)
    # print(ani)
    # print()


    # check the genome got results by ANI
    if ani:
        genus_evidence = "ani"


    # there is a ref genus assigned by ani
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

    if species and "species__" not in species:
        deepest = species
    elif genus and "genus__" not in genus:
        deepest = genus
    elif subfamily not in ["unknown", ""]:
        deepest = subfamily
    elif family not in ["unknown", ""]:
        deepest = family
    else:
        deepest = "unknown"

    # calculate the completeness proxy
    if deepest != "unknown":
        proxy = round(df.loc[genome, "length"]/taxas_lengths[deepest], 3)

        df.loc[genome, "len/taxa_len"] = proxy
        df.loc[genome, "ref_taxa"] = deepest

# -------------------------------------------------------------
# Parse taxonomy. For each reference genus, store its subfamily
lines = [line.strip().split("\t") for line in open(snakemake.params.taxonomy).readlines()[1:]]
genus_subfamily = {line[3]:line[2] for line in lines}

def assign_subfamily_sharedprots_ani(df, genus_subfamily, genome):
    genus = df.loc[genome, "genus"]
    subfamily = df.loc[genome, "subfamily"]

    if genus and "genus__" not in genus:
        if not subfamily or pd.isnull(subfamily):
            df.loc[genome, "subfamily"] = genus_subfamily[genus]






def check_no_signals_genomes(final_df, genome):

    note = str()

    # check if the genome has any assignment at family and subfamily levels
    discard = True
    check_columns = ["family", "subfamily"]
    for column in check_columns:
        value = final_df.loc[genome, column]
        if value and not pd.isnull(value):
            discard = False

    if discard:
        # remove genus & species assignments
        final_df.loc[genome, "genus"] = ""
        final_df.loc[genome, "species"] = ""
        final_df.loc[genome, "evidence_genus"] = ""


        note = "no signals"


    # check the genome has



    return discard, note





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


    # for genomes classified only by shared prots and ani, assign the subfamily
    # based on the assigned genus, if this is one of the reference genus
    assign_subfamily_sharedprots_ani(final_df, genus_subfamily, genome)


    # write completeness
    assess_completeness(final_df, taxas_lengths, genome)





    # check if there is any signal at all. Otherwise, mark the genome as discard=True
    discard, note = check_no_signals_genomes(final_df, genome)
    final_df.loc[genome, "discard"] = discard
    if note:
        notes.append(note)


    # write notes if necessary
    if notes:
        final_df.loc[genome, "notes"] = "; ".join(notes)


# Write to final output file
final_df.sort_values(by=["sample", "length", "discard"], ascending=[True, False, True], inplace=True)
final_df.to_csv(snakemake.output[0], sep="\t", header=True, index=True)
