from Bio import SearchIO

# env: utils.yaml

# Looking at hit's evalue might be too greedy, better check the ivalue of the hsp(s)

to_write = list()

for file in snakemake.input:
    records = SearchIO.parse(file, "hmmer3-text")

    for record in records:
        for hit in record.hits:
            check = False
            # if hit.is_included:
            #     to_write.append([hit.id, record.id, str(hit.evalue)])
            for hsp in hit.hsps:
                if hsp.is_included:
                    check = True

            if check:
                to_write.append([hit.id, record.id, str(hit.evalue)])


# write hits
with open(snakemake.output.hits_table, "w") as fout:
    fout.write("transeq_contig\tcrass_profile\tevalue\n")
    if to_write:
        for line in to_write:
            fout.write("\t".join(line)+ "\n")


# get the set of contigs
contigs = sorted(list(set(["_".join(line[0].split("_")[:-1]) for line in to_write])))
with open(snakemake.output.contigs, "w") as fout:
    if contigs:
        for contig in contigs:
            fout.write(f"{contig}\n")
