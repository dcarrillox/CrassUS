import pandas as pd
import glob, os
from snakemake.utils import validate
from snakemake.utils import min_version


configfile: "config/config.yaml"
report: "report/workflow.rst"


###### Parse sample sheet ######
sample_sheet = pd.read_table(config["sample_sheet"], comment='#').set_index("sample_id", drop=False)
# replace underscores in the samples_id by hyphen
sample_sheet.index = sample_sheet.index.str.replace('_','-')
#print(sample_sheet.index) # comment when running --dag

###### Wildcard constraints ######
wildcard_constraints:
    sample="|".join(sample_sheet.index),
    sample_transeq="|".join(sample_sheet.index)



###### Helper functions ######
def gather_genomes_blastall(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    crassus_fasta = expand("results/3_crass_contigs/{contig}.fasta",
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )
    ref_genomes = glob.glob("resources/genomes/*.fasta")
    return crassus_fasta + ref_genomes

def get_raw_assemblies(wildcards):
    return sample_sheet["fasta"][wildcards.sample]

def aggregate_best_codings(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    best_codings = expand(
                        "results/4_ORF/1_best_coding/{best_coding}.{ext}",
                        best_coding=glob_wildcards(f"{checkpoint_output}/{{best_coding}}.faa").best_coding,
                        ext=["faa", "gff"]
                        )
    return best_codings

def aggregate_densities(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    annotation_files = expand("results/4_ORF/0_all_codings/{contig}_{prod_ext}",
                      contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                      prod_ext=prodigal_ext
                      )
    return annotation_files

def get_prots_files(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    prots_files = expand("results/4_ORF/1_best_coding/{prots}.faa",
                      prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    return prots_files

def get_markers_files(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    markers_summary_files = expand("results/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.summary",
                      prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    markers_faa_files = expand("results/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.faa",
                      prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    markers_hmmtxt_files = expand("results/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.hmmtxt",
                      prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    return markers_summary_files + markers_faa_files + markers_hmmtxt_files

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

    tree_files = expand("results/5_phylogenies/2_trees/{marker}_trimmed.nwk",
                    marker=final_markers
                    )
    dist_files = expand("results/5_phylogenies/2_trees/{marker}_trimmed.dist",
                    marker=final_markers
                    )

    print(final_markers)
    return tree_files

def gather_dtr(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    dtr = expand("results/3_crass_contigs/dtr_blast/{contig}.dtr_blast_done",
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )
    return dtr

# the name of this function needs to be different and make clear that it aggregates
# results from several steps. "contigs" is too vague
def get_plots(wildcards):
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]

    plot = expand("results/7_ANI/2_plot/{contig}.png",
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )
    return plot

def get_genome_tables_finished(wildcards):
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    genome_tables = expand(
                        "results/4_ORF/2_functional_annot_tables/{best_coding}.table",
                        best_coding=glob_wildcards(f"{checkpoint_output}/{{best_coding}}.faa").best_coding,
                        )
    return genome_tables



def gather_gggenomes(wildcards):
    checkpoint_output = checkpoints.prepare_gggenomes_data.get(**wildcards).output[0]
    gggenomes = expand("results/7_ANI/2_plot/{gggdata}.png",
                    gggdata=glob_wildcards(f"{checkpoint_output}/{{gggdata}}.txt").gggdata,
                    )
    return gggenomes
