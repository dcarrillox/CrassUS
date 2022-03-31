#!/usr/bin/env python
# env: utils.yaml

from Bio import SearchIO

"""
Parses the results from hmmsearch to know which contigs got hits from the
Crassvirales reference set of marker genes (TerL, MCP, portal)
"""

# Link reference profile ids (such as VP02660 or VP02749) to protein names
lines = [line.strip().split("\t")
            for line in open(snakemake.params.profiles_names).readlines()]

prof_names = {line[0]:line[1] for line in lines}


# parse the results from hmmsearch. For each target (user contigs), store which
# marker profiles (queries: TerL, MCP, portal) where found.
to_write = list()
for file in snakemake.input:
    records = SearchIO.parse(file, "hmmer3-text")

    for record in records:
        for hit in record.hits:
            for hsp in hit.hsps:
                if hsp.evalue < 0.001 and hit.bitscore >= 15:
                    to_write.append([hit.id,
                                     record.id,
                                     prof_names[record.id],
                                     str(hit.evalue),
                                     str(hsp.evalue)]
                                    )
                    break



# write hits
with open(snakemake.output.hits_table, "w") as fout:
    fout.write("transeq_contig\tprofile_id\tprofile_name\te-value\ti-value\n")
    if to_write:
        for line in to_write:
            fout.write("\t".join(line)+ "\n")


# get the set of contigs
contigs = sorted(list(set(["_".join(line[0].split("_")[:-1]) for line in to_write])))
with open(snakemake.output.contigs, "w") as fout:
    if contigs:
        for contig in contigs:
            fout.write(f"{contig}\n")
