from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# env: utils.yaml

# read marker genes from file
lines = [line.strip().split("\t") for line in open("resources/yutin_2021/all_profiles/marker_profiles.txt").readlines()]
profs_names = {line[0]:line[1] for line in lines}


# parse the hmmsearch file, store hits along the profile name
markers = [marker for marker in snakemake.config["phylogenies"] if snakemake.config["phylogenies"][marker]]
markers_hits = {marker:list() for marker in markers}
records = SearchIO.parse(snakemake.input.hmmtxt, "hmmer3-text")
for record in records:
    if record.id in profs_names:
        for hit in record.hits:
            if hit.bitscore > 15:
                for hsp in hit.hsps:
                    if hsp.is_included:
                        markers_hits[profs_names[record.id]].append(hit.id)


print(markers_hits)

# read gff
gff_lines = [line.split("\t") for line in open(snakemake.input.gff).readlines() if not line.startswith("#")]
# read best coding protein file
records = {record.id:record for record in SeqIO.parse(snakemake.input.faa, "fasta")}


# iterate the names & hits
to_write_faa = list()
markers_summary = {marker:list() for marker in markers}

for marker, hits in markers_hits.items():
    # check if there were hits for that name
    # if there were hits
    if hits:
        unique_hits = list(set(hits))
        #print(marker, unique_hits)
        # init a list to store the strands to know how to merge, if necessary
        strands = list()
        # get the n_gene in the genome
        n_genes = [int(gene_id.split("|")[-1]) for gene_id in unique_hits]
        # check in the gff lines
        for n_gene in n_genes:
            # store the strand
            strands.append(gff_lines[n_gene-1][6])

        # check the presence of fragments. If so, try to merge them
        if len(unique_hits) > 1:
            # get the strand of the fragments
            strand = list(set(strands))
            if len(strand) == 1:
                # + strand, sort ascendent
                if strand[0] == "+":
                    fragments = sorted(unique_hits, key=lambda fragment: int(fragment.split("|")[-1]))
                else:
                    fragments = sorted(unique_hits, key=lambda fragment: int(fragment.split("|")[-1]), reverse=True)

                # go through the fragments. If the distance is equal or lower than 5, merge
                # the fragments. Otherwise, select the longest protein
                check = True
                for i in range(len(fragments)-1):
                    #print(fragments[i],fragments[i+1])
                    n_gene_1 = int(fragments[i].split("|")[-1])
                    n_gene_2 = int(fragments[i+1].split("|")[-1])
                    if abs(n_gene_1 - n_gene_2) > 5:
                        # get the longest protein. Sort them by their length
                        sorted_prots = sorted(fragments, key=lambda fragment: int(fragment.split("|")[-2]), reverse=True)
                        longest_prot = sorted_prots[0]
                        markers_summary[marker] = longest_prot
                        to_write_faa.append(records[longest_prot])
                        print(f"multiple copies detected: {fragments}\n{longest_prot} chosen.")
                        check = False
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
                    markers_summary[marker] = joint_id

            # if different strands, grab the longest protein
            else:
                fragments = sorted(unique_hits, key=lambda fragment: int(fragment.split("|")[-1]))
                sorted_prots = sorted(fragments, key=lambda fragment: int(fragment.split("|")[-2]), reverse=True)
                longest_prot = sorted_prots[0]
                markers_summary[marker] = longest_prot
                to_write_faa.append(records[longest_prot])
                print(f"wrong strands detected: {fragments}\n{longest_prot} chosen.")

        # otherwise, just grab the complete terminase sequence
        else:
            to_write_faa.append(records[unique_hits[0]])
            markers_summary[marker] = unique_hits[0]

    else:
        print(marker, " what")
        markers_summary[marker] = "not_found"

##
with open(snakemake.output.summary, "w") as fout:
    fout.write(f"\t" + "\t".join(markers) + "\n")
    to_write = [snakemake.wildcards.prots.split("_tbl-")[0]]
    to_write += [markers_summary[marker] for marker in markers]
    print(to_write)
    fout.write("\t".join(to_write) + "\n")

with open(snakemake.output.faa, "w") as fout:
    # if there is a TerL sequence(s) to grab from the prodigal file...
    if to_write_faa:
        SeqIO.write(to_write_faa, fout, "fasta")
