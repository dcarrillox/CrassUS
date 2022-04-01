from Bio import SeqIO
import os, glob


families_colors = {"Intestiviridae":"#EE3B3B",
                      "Crevaviridae":"#EE9A00",
                      "Suoliviridae":"#4169E1",
                      "Steigviridae":"#00CED1",
                      "Epsilon":"#CD2990",
                      "Zeta":"#006400",
                      "outgroup": "gray",
                      "NA":"gray"}


# read taxonomy ref
file = "resources/CrassUS_db/reference_taxonomy_subfamily.txt"
crass_taxonomy = dict()
lines = [line.strip().split("\t") for line in open(file).readlines()[1:]]
for line in lines:
    crass_taxonomy[line[0]] = {"family":line[1],"subfamily":line[2],"genus":line[3]}



# parse MSAs, associate a color to each genome
msas = glob.glob("resources/CrassUS_db/MSAs/*")

for msa in msas:

    to_write = ["TREE_COLORS", "SEPARATOR TAB", "DATA"]

    marker = os.path.basename(msa).replace("reference_", "").replace(".mafft-einsi", "")

    records = SeqIO.parse(msa, "fasta")

    for record in records:
        genome = record.id.split("|")[0]
        color = families_colors[crass_taxonomy[genome]["family"]]
        to_write.append(f"{record.id}\tlabel\t{color}\tbold\t2")

    with open(f"{marker}_itol.txt", "w") as fout:
        for line in to_write:
            fout.write(line + "\n")
