import os, glob
import numpy as np
from Bio import SeqIO

# env: utils.yaml

def calculate_coding_density(gff_file):
    '''
    Calculates the genome density for the given GFF
    '''
    # get contig id
    contig_id = os.path.basename(gff_file).split("_prod-")[0]
    coding = os.path.basename(gff_file).split("_prod-")[1].replace(".gff", "")

    # read GFF
    lines  = [line.strip().split("\t") for line in open(gff_file).readlines()]
    #print(lines[1])

    # get length of the genome
    length = float(lines[1][0].split(";seqlen=")[1].split(";")[0])

    # store start and ends of the ORFs
    starts = [int(line[3]) for line in lines if not line[0].startswith("#")]
    ends   = [int(line[4]) for line in lines if not line[0].startswith("#")]

    # compute average length
    lengths = list()
    for start, end in zip(starts, ends):
        lengths.append(int(end)-int(start))
    mean_len = "{:.3f}".format(np.mean(lengths))

    # compute coding length. It collapses overlaping ORFs in the same
    # interval
    intervals = [[s,e] for s,e in zip(starts, ends)]

    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)

    # coding length
    coding_bases = float(0)
    for interval in merged:
        coding_bases += 1 + (interval[1] - interval[0])

    # coding density
    density = "{:.3f}".format(coding_bases/length)

    # return id, 11/TAG/TGA, seqlen, density, n_prots, prots_len
    return [contig_id, coding, str(int(length)), density, str(len(starts)), mean_len]

def main():

    gffs_files = [file for file in snakemake.input if file.endswith(".gff")]

    codings = list()
    for file in gffs_files:
        to_add = calculate_coding_density(file)
        codings.append(to_add)


    # pick the best contig
    contigs_codings = {line[0]:list() for line in codings}
    for line in codings:
        contigs_codings[line[0]].append(line)

    to_write = list()
    for contig, codings in contigs_codings.items():
        sorted_codings = sorted(codings, key=lambda coding: float(coding[3]), reverse=True)

        # check density & mean length, pick the best one
        dens1, dens2 = float(sorted_codings[0][3]), float(sorted_codings[1][3])
        leng1, leng2 = float(sorted_codings[0][5]), float(sorted_codings[1][5])
        if dens1 - dens2 > 0.05 or leng1 - leng2 > 200:
            sorted_codings[0].append("<--")
        else:
            for coding in sorted_codings:
                if coding[1] == "11":
                    coding.append("<--")
        # sort by coding name
        sorted_codings = sorted(sorted_codings, key=lambda coding: coding[1])
        for coding in sorted_codings:
            to_write.append(coding)

    with open(snakemake.output[0], "w") as fout:
        fout.write("contig\tcoding\tlength\tdensity\tn_prots\tmean_len\n")
        for line in to_write:
            fout.write("\t".join(line) + "\n")




if __name__ == "__main__":
    main()
