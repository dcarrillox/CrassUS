from Bio import SeqIO

records = SeqIO.parse(snakemake.input[0], "fasta")

# measure the length of each record, store in a dict
ids_records = {record.id:record for record in records}

# sort by value (length)
sorted_ids_records = {id: record for id, record in sorted(ids_records.items(), key=lambda item: len(item[1].seq), reverse=True)}

# replace the header by {sample_contigX_length}
to_write_fasta = list()
to_write_table = list()

n = 1
for id, record in sorted_ids_records.items():
    to_add_table = [record.id]
    if len(record.seq) >= 5000:
        rehead = f"{snakemake.wildcards.sample}_{n}_{len(record.seq)}"
        to_add_table.append(rehead)
        record.id = rehead
        record.description = ""
        to_write_fasta.append(record)
        n += 1

    else:
        to_add_table.append("")

    to_write_table.append(to_add_table)


with open(snakemake.output.fasta, "w") as fout:
    SeqIO.write(to_write_fasta, fout, "fasta")

with open(snakemake.output.table, "w") as fout:
    fout.write("old_id\tnew_id\n")
    for line in to_write_table:
        fout.write("\t".join(line) + "\n")
