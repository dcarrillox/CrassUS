from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# Looking at hit's evalue might be too greedy, better check the ivalue of the hsp(s)
records = SearchIO.parse(snakemake.input.hmmtxt, "hmmer3-text")

terL_profiles = ["VP00412", "VP02579", "VP02686"]

# identify hits to the terminase profiles
terL_hits = list()
for record in records:
    if record.id in terL_profiles:
        for hit in record.hits:
            for hsp in hit.hsps:
                if hsp.is_included:
                    terL_hits.append(hit.id)

# check if any terL hits were identified
to_write_summary = list()

if terL_hits:
    # make the hits uniq
    terL_hits = list(set(terL_hits))
    # init a list to store the ids for the faa sequence
    terL_faa  = list()

    # read gff file
    gff_lines = [line.split("\t") for line in open(snakemake.input.gff).readlines() if not line.startswith("#")]

    # check that line gene_id - 1 in the gff is not partial
    to_write_summary_summary = list()
    for gene in terL_hits:
        n_gene = int(gene.split("|")[-1])
        if "partial=00" in gff_lines[n_gene-1][-1]:
            to_write_summary.append(f"{gene}\t{gff_lines[n_gene-1][6]}\tTerL")
            # gather gene_id, n_gene and strand, I will need in case of merging terL fragments
            terL_faa.append([gene, int(gene.split("|")[-1]), gff_lines[n_gene-1][6]])
        else:
            to_write_summary.append(f"{gene}\t{gff_lines[n_gene-1][6]}\tcontig end TerL")

else:
    to_write_summary("No TerL hits")


# write summary file
with open(snakemake.output.summary, "w") as fout:
    for line in to_write_summary:
        fout.write(f"{line}\n")


# check if there are .faa proteins, and if they have to be merged in some way
to_write_faa = list()
if terL_faa:
    # read prodigal .faa file
    records = {record.id:record for record in SeqIO.parse(snakemake.input.faa, "fasta")}

    # check the presence of fragments. If so, try to merge them
    if len(terL_faa) > 1:
        # get the strand of the terL fragments
        strand = list(set([gene[2] for gene in terL_faa]))
        if len(strand) == 1:
            # + strand, sort ascendent
            if strand[0] == "+":
                terL_fragments = sorted(terL_faa, key=lambda fragment: fragment[1])
            else:
                terL_fragments = sorted(terL_faa, key=lambda fragment: fragment[1], reverse=True)


        # go through the fragments. If the distance is equal or lower than 3, merge
        # the fragments. Otherwise, abort it.
        check = True
        print(strand)
        print(terL_fragments)
        for i in range(len(terL_fragments)-1):
            print(terL_fragments[i],terL_fragments[i+1])
            if abs(terL_fragments[i][1] - terL_fragments[i+1][1]) > 3:
                check = False
                break

        if check:
            contig = terL_fragments[0][0].split("|")[0]
            joint_n = "_".join([str(fragment[1]) for fragment in terL_fragments])
            seq = str()
            for fragment in terL_fragments:
                seq += str(records[fragment[0]].seq).replace("*", "")

            joint_id = f"{contig}|{len(seq)}|{joint_n}"
            new_record = SeqRecord(Seq(seq), id=joint_id, description="")

            to_write_faa.append(new_record)


    # otherwise, just grab the complete terminase sequence
    else:
        to_write_faa.append(records[terL_faa[0][0]])






# write .faa file
with open(snakemake.output.faa, "w") as fout:
    # if there is a TerL sequence(s) to grab from the prodigal file...
    if to_write_faa:
        SeqIO.write(to_write_faa, fout, "fasta")
