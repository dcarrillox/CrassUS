import pandas as pd
import glob, os
from snakemake.utils import validate
from snakemake.utils import min_version


configfile: "config/config.yaml"
report: "report/workflow.rst"


###### Parse sample sheet ######
sample_sheet = pd.read_table(config["sample_sheet"]).set_index("sample_id", drop=False)

###### Wildcard constraints ######
wildcard_constraints:
    sample="|".join(sample_sheet.index),
    sample_transeq="|".join(sample_sheet.index)



###### Helper functions ######

def get_raw_assemblies(wildcards):
    return sample_sheet["fasta"][wildcards.sample]

def get_20k_assemblies(wildcards):
    return checkpoints.filter_assemblies.get(sample=wildcards.sample).output

def aggregate_best_codings(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    best_codings = expand(
                        "results/4_prodigal/best_coding/{best_coding}.{ext}",
                        best_coding=glob_wildcards(f"{checkpoint_output}/{{best_coding}}.faa").best_coding,
                        ext=["faa", "gff"]
                        )
    return best_codings

def aggregate_pyani(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]

def aggregate_densities(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    prodigal_files = expand("results/4_prodigal/all_codings/{contig}_{prod_ext}",
                      contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                      prod_ext=prodigal_ext
                      )
    return prodigal_files

def get_prots_files(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    prots_files = expand("results/4_prodigal/best_coding/{prots}.faa",
                      prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    return prots_files

def get_markers_files(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    markers_summary_files = expand("results/5_phylogenies/markers_scan/{prots}_markers.summary",
                      prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    markers_faa_files = expand("results/5_phylogenies/markers_scan/{prots}_markers.faa",
                      prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    return markers_summary_files + markers_faa_files

def gather_trees(wildcards):
    # check the config file to know which markers are requested for the trees
    config_markers = list()
    for marker in config["phylogenies"]:
        if config["phylogenies"][marker]:
            config_markers.append(marker)

    # get the markers given by the checkpoint
    checkpoint_output = checkpoints.summarize_markers.get(**wildcards).output.faa_dir
    checkpoint_markers = glob_wildcards(f"{checkpoint_output}/{{marker}}.faa").marker

    # as final, take the intersection between the given by the checpoint and the requested from the config file
    final_markers = [marker for marker in checkpoint_markers if marker in config_markers]
    mafft_files = expand("results/5_phylogenies/msa/{marker}.msa",
                    marker=final_markers
                    )

    print(final_markers)
    return mafft_files


# the name of this function needs to be different and make clear that it aggregates
# results from several steps. "contigs" is too vague
def aggregate_contigs(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    blast    = expand("results/5_blast/{contig}.blast",
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )
    prodigal = expand("results/4_prodigal/all_codings/{contig}_{prod_ext}",
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
    #return prodigal + fastani
    return prodigal
