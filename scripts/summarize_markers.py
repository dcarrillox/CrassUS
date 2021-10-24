import pandas as pd
from Bio import SeqIO, SearchIO
import os

# env: utils.yaml


# parse all the hmmscan results and save to df
hmmtxt_files = [file for file in snakemake.input if file.endswith(".hmmtxt")]
columns = ["TerL_prot", "TerL_eval", "TerL_cov", "TerL_prot_cov",
           "MCP_prot", "MCP_eval", "MCP_cov", "MCP_prot_cov",
           "portal_prot", "portal_eval", "portal_cov", "portal_prot_cov"]
hmmscan_df = pd.DataFrame(columns=columns)
# profiles length for their cov calculation. Get also the type
profiles_length = {line.split("\t")[0]:int(line.strip().split("\t")[2]) for line in open(snakemake.params.profiles_length).readlines()}
profiles_type = {line.split("\t")[0]:line.strip().split("\t")[1] for line in open(snakemake.params.profiles_length).readlines()}


for hmmtxt_file in hmmtxt_files:
    # check it is not an empty file
    if os.path.getsize(hmmtxt_file) != 0:
        records = SearchIO.parse(hmmtxt_file, "hmmer3-text")

        for record in records:
            contig = record.id.split("|")[0]
            protein_starts = list()
            protein_ends   = list()
            hit_starts = list()
            hit_ends   = list()
            # break after the first hit, they are sorted by evalue
            for hit in record.hits:
                hit_len = profiles_length[hit.id] # run cell below to get this dictionary. hmmsearch is the only way I have to retrieve the length of the hmm
                marker = profiles_type[hit.id]
                hit_evalue =  hit.evalue

                for hsp in hit.hsps:
                    hit_starts.append(hsp.query_start)
                    hit_ends.append(hsp.query_end)
                    protein_starts.append(hsp.env_start)
                    protein_ends.append(hsp.env_end)
                break

            # protein coverage
            prot_len = int(record.id.split("|")[1])
            intervals = [[s,e] for s,e in zip(protein_starts, protein_ends)]
            intervals.sort(key=lambda interval: interval[0])
            #print(record.id, "\n",hit)
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

            prot_cov = round(float(covered_length/prot_len),3)


            # marker coverage
            intervals = [[s,e] for s,e in zip(hit_starts, hit_ends)]

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

            marker_cov = round(float(covered_length/hit_len),3)

            # fill the hmmscan_df in
            hmmscan_df.loc[contig, f"{marker}_prot"] = record.id
            hmmscan_df.loc[contig, f"{marker}_eval"] = hit_evalue
            hmmscan_df.loc[contig, f"{marker}_cov"] = marker_cov
            hmmscan_df.loc[contig, f"{marker}_prot_cov"] = prot_cov



# gather all the summary files and put them in a single table
summary_files = [file for file in snakemake.input if file.endswith(".summary")]
rows = list()
for summary_file in summary_files:
    lines = [line.strip().split("\t") for line in open(summary_file).readlines()]
    rows.append(lines[-1])

columns = ["contig", "TerL", "MCP", "portal"]
summary_df = pd.DataFrame(rows, columns=columns)
summary_df = summary_df.sort_values(by="contig")
summary_df = summary_df.set_index("contig")


# store all the .faa files in the all_records dictionary
faa_files = [file for file in snakemake.input if file.endswith(".faa") and os.path.getsize(file) > 0]
all_records = dict()
for faa_file in faa_files:
    records = SeqIO.parse(faa_file, "fasta")
    for record in records:
        all_records[record.id] = record


# create the output dir for the .faa files
os.makedirs(snakemake.output.faa_dir, exist_ok = False)

# iterate the contigs while checking their markers in the hmmscan_df, and write to final marker.faa file
for marker in summary_df.columns: # don't look at the "contigs" column
    #print(f"marker -> {marker}")
    to_write = list()
    for contig in summary_df.index:
        #print(f"\tcontig -> {contig}")
        # process only contigs for which there were proteins passing the multiple copies etc conditions, for this specific marker.
        if contig in hmmscan_df.index:
            if not pd.isnull(hmmscan_df.loc[contig, f"{marker}_prot"]):
                print(hmmscan_df.loc[contig, f"{marker}_prot"])
                prot_id = hmmscan_df.loc[contig, f"{marker}_prot"]
                marker_cov = hmmscan_df.loc[contig, f"{marker}_cov"]
                prot_cov   = hmmscan_df.loc[contig, f"{marker}_prot_cov"]
                if marker_cov >= 0.5 and prot_cov >= 0.5:
                    to_write.append(all_records[prot_id])
                else:
                    #print(prot_id, marker_cov, prot_cov, hmmscan_df.loc[contig, f"{marker}_eval"])
                    summary_df.loc[contig, marker] = "too_short"

    #print(to_write)
    if to_write:
        with open(f"{snakemake.output.faa_dir}/{marker}.faa", "w") as fout:
            SeqIO.write(to_write, fout, "fasta")


# write the final markers.summary file and the coverages file
summary_df.to_csv(snakemake.output.summary, sep="\t")
hmmscan_df.sort_index(inplace=True)
hmmscan_df.to_csv(snakemake.output.coverages, sep="\t")
