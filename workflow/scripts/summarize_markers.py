import pandas as pd
from Bio import SeqIO, SearchIO
import os

# env: utils.yaml


def calculate_coverage(length, starts, ends):
    intervals = [[s,e] for s,e in zip(starts, ends)]
    intervals.sort(key=lambda interval: interval[0])

    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    # calculate coverage
    covered_length = float(0)
    for interval in merged:
        covered_length += 1 + (interval[1] - interval[0])

    cov = round(float(covered_length/length),3)

    return cov


# -------------------------------
# Set a DataFrame with columns "_prot", "_eval", "_cov" and "_prot_cov" for each
# marker.
markers = sorted([marker for marker in snakemake.config["phylogenies"]["markers"] if snakemake.config["phylogenies"]["markers"][marker]], reverse=True)
columns = ["prot", "eval", "cov", "prot_cov"]
final_columns = list()
for marker in markers:
    for column in columns:
        final_columns.append(f"{marker}_{column}")

hmmscan_df = pd.DataFrame(columns=final_columns)



# ------------------------------------------------------------------
# Parse hmmscan files for all the contigs while calculating coverage. Fill the
# df created above.

# get name and length of the profiles
lines = [line.strip().split("\t") for line in open(snakemake.params.profiles_length).readlines()]
profiles_length = {line[0]:int(line[2]) for line in lines}
profiles_type = {line[0]:line[1] for line in lines}

# iterate hmmscan files
hmmtxt_files = [file for file in snakemake.input if file.endswith(".hmmtxt")]
for hmmtxt_file in hmmtxt_files:
    # check it is not an empty file
    if os.path.getsize(hmmtxt_file) != 0:
        records = SearchIO.parse(hmmtxt_file, "hmmer3-text")

        # store metrics to calculate coverage
        for record in records:
            contig = record.id.split("|")[0]
            protein_starts = list()
            protein_ends   = list()
            marker_starts = list()
            marker_ends   = list()
            # break after the first hit, they are sorted by evalue
            for hit in record.hits:
                marker_length = profiles_length[hit.id]
                marker = profiles_type[hit.id]
                hit_evalue =  hit.evalue

                for hsp in hit.hsps:
                    marker_starts.append(hsp.query_start)
                    marker_ends.append(hsp.query_end)
                    protein_starts.append(hsp.env_start)
                    protein_ends.append(hsp.env_end)
                break


            # Calculate protein coverage
            protein_length = int(record.id.split("|")[1])
            prot_cov = calculate_coverage(protein_length, protein_starts, protein_ends)

            # Calculate marker coverage
            marker_cov = calculate_coverage(marker_length, marker_starts, marker_ends)


            # fill the hmmscan_df in
            hmmscan_df.loc[contig, f"{marker}_prot"] = record.id
            hmmscan_df.loc[contig, f"{marker}_eval"] = hit_evalue
            hmmscan_df.loc[contig, f"{marker}_cov"] = marker_cov
            hmmscan_df.loc[contig, f"{marker}_prot_cov"] = prot_cov



# --------------------------------------------------------
# store all the markers FAA files in the "all_records" dictionary
faa_files = [file for file in snakemake.input if file.endswith(".faa") and os.path.getsize(file) > 0]
all_records = dict()
for faa_file in faa_files:
    records = SeqIO.parse(faa_file, "fasta")
    for record in records:
        all_records[record.id] = record



# ------------------------------------------------------------------------------
# gather all the summary files and put them in a single table "summary_df". Used
# later to check that the marker was found in the contig.
summary_files = [file for file in snakemake.input if file.endswith(".summary")]
rows = list()
for summary_file in summary_files:
    lines = [line.strip().split("\t") for line in open(summary_file).readlines()]
    rows.append(lines[-1])

columns = ["contig"] + markers
summary_df = pd.DataFrame(rows, columns=columns)
summary_df = summary_df.sort_values(by="contig")
summary_df = summary_df.set_index("contig")



# ------------------------------------------------------------------------------
# iterate the markers found on each contig checking their coverage. If they pass
# the cutoff, write to final FAA file. Otherwise,

os.makedirs(snakemake.output.faa_dir, exist_ok = True)

for marker in summary_df.columns: # don't look at the "contigs" column
    to_write = list()
    for contig in summary_df.index:
        # process only contigs for which there were proteins passing the multiple
        # copies etc conditions, for this specific marker.
        if contig in hmmscan_df.index:
            if not pd.isnull(hmmscan_df.loc[contig, f"{marker}_prot"]):
                prot_id = hmmscan_df.loc[contig, f"{marker}_prot"]
                marker_cov = hmmscan_df.loc[contig, f"{marker}_cov"]
                prot_cov   = hmmscan_df.loc[contig, f"{marker}_prot_cov"]

                cutoff_marker_cov = snakemake.config["phylogenies"]["marker_cov"]
                cutoff_protein_cov = snakemake.config["phylogenies"]["marker_cov"]

                if marker_cov >= cutoff_marker_cov and prot_cov >= cutoff_protein_cov:
                    to_write.append(all_records[prot_id])
                else:
                    summary_df.loc[contig, marker] = "too_short"

    # write final FAA file of the marker
    if to_write:
        with open(f"{snakemake.output.faa_dir}/{marker}.faa", "w") as fout:
            SeqIO.write(to_write, fout, "fasta")



# -----------------------------------------------------------
# write the final markers.summary file and the coverages file
summary_df.to_csv(snakemake.output.summary, sep="\t")
hmmscan_df.sort_index(inplace=True)
hmmscan_df.to_csv(snakemake.output.coverages, sep="\t")
