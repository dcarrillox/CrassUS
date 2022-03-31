#!/usr/bin/env python
# env: utils.yaml

from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys


# -----------------------------------------------------------------
# Find out which contig genes (hit.id) were annotated as any of the markers

# store marker_name - marker_profile_id pairs
lines = [line.strip().split("\t") for line in open(snakemake.params.markers_ids).readlines()]
profs_names = {line[0]:line[1] for line in lines}

# get from the config file which markers were used, init a dict with them
markers = sorted([marker for marker in snakemake.config["phylogenies"]["markers"] if snakemake.config["phylogenies"]["markers"][marker]], reverse=True)
markers_hits = {marker:list() for marker in markers}

# parse the hmmsearch file, store hits along the marker name
records = SearchIO.parse(snakemake.input.hmmtxt, "hmmer3-text")
for record in records:
    if record.id in profs_names:
        for hit in record.hits:
            for hsp in hit.hsps:
                if hsp.evalue < 0.001 and hit.bitscore >= 15:
                    markers_hits[profs_names[record.id]].append(hit.id)



# ---------------------------------------------------------------------
# For the best coding, read GFF (strand info) and FAA (sequence itself)
gff_lines = [line.split("\t") for line in open(snakemake.input.gff).readlines() if not line.startswith("#")]
records = {record.id:record for record in SeqIO.parse(snakemake.input.faa, "fasta")}



# --------------------------------------------------------
# iterate the markers and their hits, merging if necessary
to_write_faa = list()
markers_summary = {marker:list() for marker in markers}

for marker, hits in markers_hits.items():
    # check if there were hits for that name
    if not hits:
        markers_summary[marker] = "not_found"
    else:
        # remove redundancy: hits with many profiles of the marker
        unique_hits = list(set(hits))
        # init a list to store the strands to know how to merge, if necessary
        strands = list()
        # get the n_gene in the genome
        n_genes = [int(gene_id.split("|")[-1]) for gene_id in unique_hits]
        # check in the gff lines
        for n_gene in n_genes:
            # store the strand
            strands.append(gff_lines[n_gene-1][6])

        # check the presence of fragments. If so, try to merge them
        if len(unique_hits) == 1: # only one fragment
            to_write_faa.append(records[unique_hits[0]])
            markers_summary[marker] = unique_hits[0]
        else: # multiple fragments
            # get the strand of the fragments
            strand = list(set(strands))
            if len(strand) != 1: # fragments on different strands
                fragments = sorted(unique_hits, key=lambda fragment: int(fragment.split("|")[-1]))
                sorted_prots = sorted(fragments, key=lambda fragment: int(fragment.split("|")[-2]), reverse=True)
                longest_prot = sorted_prots[0]
                markers_summary[marker] = longest_prot
                to_write_faa.append(records[longest_prot])
                print(f"wrong strands detected: {fragments}\n{longest_prot} chosen.",file=sys.stderr)

            else: # fragments on the same strand
                # + strand, sort ascendent
                if strand[0] == "+":
                    fragments = sorted(unique_hits, key=lambda fragment: int(fragment.split("|")[-1]))
                else:
                    fragments = sorted(unique_hits, key=lambda fragment: int(fragment.split("|")[-1]), reverse=True)

                # go through the fragments. If the distance is <= 5 genes, merge
                # the fragments. Otherwise, select the longest one
                check = True
                for i in range(len(fragments)-1):
                    n_gene_1 = int(fragments[i].split("|")[-1])
                    n_gene_2 = int(fragments[i+1].split("|")[-1])
                    if abs(n_gene_1 - n_gene_2) > 5:
                        # get the longest protein. Sort them by their length
                        sorted_prots = sorted(fragments, key=lambda fragment: int(fragment.split("|")[-2]), reverse=True)
                        longest_prot = sorted_prots[0]
                        markers_summary[marker] = longest_prot
                        to_write_faa.append(records[longest_prot])
                        print(f"multiple copies detected: {fragments}\n{longest_prot} chosen.", file=sys.stderr)
                        check = False
                        break

                if check:
                    contig = fragments[0].split("|")[0]
                    combi_n = "_".join([fragment.split("|")[-1] for fragment in fragments])
                    seq = str()
                    for fragment in fragments:
                        seq += str(records[fragment].seq).replace("*", "")

                    combi_id = f"{contig}|{len(seq)}|{combi_n}"
                    new_record = SeqRecord(Seq(seq), id=combi_id, description="")

                    to_write_faa.append(new_record)
                    markers_summary[marker] = combi_id




# ----------------------
# Write final FAA file and summary

with open(snakemake.output.faa, "w") as fout:
    if to_write_faa:
        SeqIO.write(to_write_faa, fout, "fasta")

with open(snakemake.output.summary, "w") as fout:
    fout.write(f"\t" + "\t".join(markers) + "\n")
    to_write = [snakemake.wildcards.prots.split("_tbl-")[0]]
    to_write += [markers_summary[marker] for marker in markers]
    fout.write("\t".join(to_write) + "\n")
