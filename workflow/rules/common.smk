import pandas as pd
import glob, os
from snakemake.utils import validate
from snakemake.utils import min_version


configfile: "config/config.yaml"
report: "report/workflow.rst"


###### Parse sample sheet ######
sample_sheet = pd.read_table(config["samples"]).set_index("sample_id", drop=False)

# check if there are fastq samples in the input
in_fastq = False
if "fastq" in sample_sheet["type"].values:
    in_fastq = True

# check if there are assemblies in the input
in_assembly = False
if "assembly" in sample_sheet["type"].values:
    in_assembly = True



###### Wildcard constraints ######
wildcard_constraints:
    sample="|".join(sample_sheet.index),
    sample_transeq="|".join(sample_sheet.index)



###### Helper functions ######

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = sample_sheet.loc[(wildcards.sample), ["file1", "file2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.file1, "r2": fastqs.file2}
    return {"r1": fastqs.file1}

def is_single_end(sample):
    """Return True if sample is single end."""
    return pd.isnull(sample_sheet.loc[(sample), "file2"])

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample."""
    if not is_single_end(wildcards.sample):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}_{group}.fastq",
            sample=wildcards.sample,
            group=["R1", "R2"])

    return [f"results/trimmed/{wildcards.sample}.fastq"]

def get_final_reads(wildcards):
    # pre-processing True
    if config["fastq_processing"]:
        if not is_single_end(wildcards.sample):
            a= expand("results/trimmed/{sample}_{group}.fastq",
            sample=wildcards.sample, group=["R1", "R2"])
            return a

        b = [f"results/trimmed/{wildcards.sample}.fastq"]
        return b

    # pre-processing False
    fastqs = sample_sheet.loc[(wildcards.sample), ["file1", "file2"]].dropna().tolist()
    return fastqs

def get_raw_assemblies(wildcards):
    if sample_sheet["type"][wildcards.sample] == "fastq":
        return "results/1_assembly/0_raw/{sample}/contigs.fasta"
    else:
        return sample_sheet["file1"][wildcards.sample]

def get_20k_assemblies(wildcards):
    return checkpoints.filter_assemblies.get(sample=wildcards.sample).output

def format_spades_input(wildcards):
    # pre-processing True
    if config["fastq_processing"]:
        if not is_single_end(wildcards.sample):
            a= f"--meta -1 results/trimmed/{wildcards.sample}_R1.fastq -2 results/trimmed/{wildcards.sample}_R2.fastq"
            return a

        b = f"-s results/trimmed/{wildcards.sample}.fastq"
        return b

    # pre-processing False
    if not is_single_end(wildcards.sample):
        a= f"--meta -1 {sample_sheet['file1'][wildcards.sample]} -2 {sample_sheet['file2'][wildcards.sample]}"
        return a

    b = f"-s {sample_sheet['file1'][wildcards.sample]}"
    return b

def aggregate_best_codings(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    best_codings = expand(
                        "results/5_prodigal/best_coding/{best_coding}.{ext}",
                        best_coding=glob_wildcards(f"{checkpoint_output}/{{best_coding}}.faa").best_coding,
                        ext=["faa", "gff"]
                        )
    return best_codings

def aggregate_pyani(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]

def aggregate_densities(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    prodigal_files = expand("results/5_prodigal/all_codings/{contig}_{prod_ext}",
                      contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                      prod_ext=prodigal_ext
                      )
    return prodigal_files

def get_prots_files(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    prots_files = expand("results/5_prodigal/best_coding/{prots}.faa",
                      prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    return prots_files

def aggregate_contigs(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    blast    = expand("results/5_blast/{contig}.blast",
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )
    prodigal = expand("results/5_prodigal/all_codings/{contig}_{prod_ext}",
                      contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                      prod_ext=prodigal_ext
                      )
    fastani =  expand("results/6_fastani/{contig}.fastani",
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )

    pyani = expand("results/7_pyani/{contig}/.done",
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )

    #return prodigal + pyani
    return prodigal + fastani
