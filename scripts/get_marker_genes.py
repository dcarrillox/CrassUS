from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# read marker genes from file
lines = [line.strip().split("\t") for line in open("resources/yutin_2021/all_profiles/marker_profiles.txt").readlines()]
profs_names = {line[0]:line[1] for line in lines}

# parse the hmmsearch file, store hits along the profile name
names = ["TerL", "MCP", "portal", "primase"]
names_hits = {name:list() for name in names}
records = SearchIO.parse(snakemake.input.hmmtxt, "hmmer3-text")
for record in records:
    if record.id in profs_names:
        for hit in record.hits:
            if hit.bitscore > 15:
                for hsp in hit.hsps:
                    if hsp.is_included:
                        names_hits[profs_names[record.id]].append(hit.id)


#print(names_hits)

# read gff
gff_lines = [line.split("\t") for line in open(snakemake.input.gff).readlines() if not line.startswith("#")]
# read best coding protein file
records = {record.id:record for record in SeqIO.parse(snakemake.input.faa, "fasta")}


# iterate the names & hits
to_write_faa = list()
names_summary = {name:list() for name in names}

for name, hits in names_hits.items():
    # check if there were hits for that name
    # if there were hits
    if hits:
        hits = list(set(hits))
        # init a list to store the strands to know how to merge, if necessary
        strands = list()
        # check if any of them is truncated in the edges of the contig
        truncated = False
        # get the n_gene in the genome
        n_genes = [int(gene_id.split("|")[-1]) for gene_id in hits]
        # check in the gff lines
        for n_gene in n_genes:
            # store the strand
            strands.append(gff_lines[n_gene-1][6])
            if "partial=00" not in gff_lines[n_gene-1][-1]:
                truncated = True

        if not truncated:
            # check the presence of fragments. If so, try to merge them
            if len(hits) > 1:
                # get the strand of the fragments
                strand = list(set(strands))
                if len(strand) == 1:
                    # + strand, sort ascendent
                    if strand[0] == "+":
                        fragments = sorted(hits, key=lambda fragment: int(fragment.split("|")[-1]))
                    else:
                        fragments = sorted(hits, key=lambda fragment: int(fragment.split("|")[-1]), reverse=True)


                # go through the fragments. If the distance is equal or lower than 3, merge
                # the fragments. Otherwise, abort it.
                check = True
                print(strand, snakemake.wildcards.prots, name )
                print(fragments)
                for i in range(len(fragments)-1):
                    print(fragments[i],fragments[i+1])
                    n_gene_1 = int(fragments[i].split("|")[-1])
                    n_gene_2 = int(fragments[i+1].split("|")[-1])
                    if abs(n_gene_1 - n_gene_2) > 5:
                        check = False
                        names_summary[name] = "multiple copies"
                        break

                if check:
                    contig = fragments[0].split("|")[0]
                    joint_n = "_".join([fragment.split("|")[-1] for fragment in fragments])
                    seq = str()
                    for fragment in fragments:
                        seq += str(records[fragment].seq).replace("*", "")

                    joint_id = f"{contig}|{len(seq)}|{joint_n}"
                    new_record = SeqRecord(Seq(seq), id=joint_id, description="")

                    to_write_faa.append(new_record)
                    names_summary[name] = joint_id


            # otherwise, just grab the complete terminase sequence
            else:
                to_write_faa.append(records[hits[0]])
                names_summary[name] = hits[0]


        else:
            names_summary[name] = "truncated"

    else:
        names_summary[name] = "Not found"

print(names_summary)
##
with open(snakemake.output.summary, "w") as fout:
    fout.write(f"\t" + "\t".join(names) + "\n")
    to_write = [snakemake.wildcards.prots.split("_prod")[0]]
    to_write += [names_summary[name] for name in names]
    fout.write("\t".join(to_write) + "\n")

with open(snakemake.output.faa, "w") as fout:
    # if there is a TerL sequence(s) to grab from the prodigal file...
    if to_write_faa:
        SeqIO.write(to_write_faa, fout, "fasta")
