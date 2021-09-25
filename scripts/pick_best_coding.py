from Bio import SeqIO
import glob,os

# env: utils.yaml

def format_faa_headers(faa_file, out_dir):
    contig_id = os.path.basename(faa_file).split("_tbl-")[0]

    records = SeqIO.parse(faa_file, "fasta")

    to_write= list()
    for record in records:
        record.description = ""
        nprot = record.id.split("_")[-1]
        record.id = f"{contig_id}|{len(record.seq)}|{nprot}"
        to_write.append(record)

    outfile_path = os.path.join(out_dir, os.path.basename(faa_file))
    with open(outfile_path, "w") as fout:
        SeqIO.write(to_write, fout, "fasta")



def main():


    os.makedirs(snakemake.output[0], exist_ok=True)

    # copy .faa files of best codings to the root folder
    best_codings = [line.strip().split("\t") for line in open(snakemake.input[0]).readlines() if line.endswith("<--\n")]
    for coding in best_codings:
        raw_prefix = f"{snakemake.params.raw_dir}/{coding[0]}_tbl-{coding[1]}"

        # copy the .gff file
        os.system(f"cp {raw_prefix}.gff {snakemake.output[0]}")

        # read faa file and remove the header description
        in_faa_file = f"{raw_prefix}.faa"
        format_faa_headers(in_faa_file, snakemake.output[0])


if __name__ == "__main__":
    main()
