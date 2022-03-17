import pandas as pd
import glob, os, sys
from snakemake.utils import validate
from snakemake.utils import min_version


configfile: "config/config.yaml"
report: "report/workflow.rst"


###### Parse sample sheet ######
sample_sheet = pd.read_table(config["sample_sheet"], comment='#').set_index("sample_id", drop=False)
# replace underscores in the samples_id by hyphen
sample_sheet.index = sample_sheet.index.str.replace('_','-')
# check there is only one analysis_id
ANALYSES_IDS = list(set(sample_sheet.analysis_id.tolist()))
if len(ANALYSES_IDS) > 1:
    print(f"More than one 'analysis_id' in the sample sheet '{os.path.basename(config['sample_sheet'])}':")
    for id in ANALYSES_IDS:
        print(f"\t- {id}")
    print("Please run only one analysis each time. Exiting...")
    sys.exit()
else:
    SAMPLES=sample_sheet.index.tolist()
#print(sample_sheet.index) # comment when running --dag


###### Wildcard constraints ######
wildcard_constraints:
    analysis_id="|".join(ANALYSES_IDS),
    sample="|".join(sample_sheet.index),
    sample_transeq="|".join(sample_sheet.index),


########## Check reference db ############
if not os.path.isfile("./resources/crassus_dependencies/version"):
    print("Setting up reference database and dependencies...")

    # replace by Zenodo download
    os.system("cp -r crassus_dependencies/ resources")
    print("Done")



###### Helper functions ######
def get_raw_assemblies(wildcards): # used
    return sample_sheet["fasta"][wildcards.sample]

def aggregate_densities(wildcards): #used
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    annotation_files = expand("results/{analysis_id}/4_ORF/0_all_codings/{contig}_{prod_ext}",
                              analysis_id=ANALYSES_IDS[0],
                              contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                              prod_ext=prodigal_ext
                              )
    return annotation_files

def get_prots_files(wildcards): # used
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    prots_files = expand("results/{analysis_id}/4_ORF/1_best_coding/{prots}.faa",
                        analysis_id=ANALYSES_IDS[0],
                        prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                      )
    return prots_files

def get_genome_tables_finished(wildcards): # used
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    genome_tables = expand(
                        "results/{analysis_id}/4_ORF/2_functional_annot_tables/{best_coding}.table",
                        analysis_id=ANALYSES_IDS[0],
                        best_coding=glob_wildcards(f"{checkpoint_output}/{{best_coding}}.faa").best_coding,
                        )
    return genome_tables

def get_markers_files(wildcards): # used
    checkpoint_output = checkpoints.pick_best_coding.get(**wildcards).output[0]
    markers_summary_files = expand("results/{analysis_id}/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.summary",
                            analysis_id=ANALYSES_IDS[0],
                            prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                            )
    markers_faa_files = expand("results/{analysis_id}/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.faa",
                        analysis_id=ANALYSES_IDS[0],
                        prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                        )
    markers_hmmtxt_files = expand("results/{analysis_id}/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.hmmtxt",
                            analysis_id=ANALYSES_IDS[0],
                            prots=glob_wildcards(f"{checkpoint_output}/{{prots}}.faa").prots
                            )
    return markers_summary_files + markers_faa_files + markers_hmmtxt_files

def gather_genomes_blastall(wildcards): # used
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    crassus_fasta = expand("results/{analysis_id}/3_crass_contigs/{contig}.fasta",
                    analysis_id=ANALYSES_IDS[0],
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )
    ref_genomes = glob.glob("resources/crassus_dependencies/reference_genomes/sequences/*.fasta")
    return crassus_fasta + ref_genomes

def gather_trees(wildcards): # used
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

    tree_files = expand("results/{analysis_id}/5_phylogenies/2_trees/{marker}_trimmed.nwk",
                    analysis_id=ANALYSES_IDS[0],
                    marker=final_markers
                    )
    dist_files = expand("results/{analysis_id}/5_phylogenies/2_trees/{marker}_trimmed.dist",
                    analysis_id=ANALYSES_IDS[0],
                    marker=final_markers
                    )

    print(final_markers)
    return tree_files

def gather_dtr(wildcards): # used
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    dtr = expand("results/{analysis_id}/3_crass_contigs/dtr_blast/{contig}.dtr_blast_done",
                    analysis_id=ANALYSES_IDS[0],
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )
    return dtr

def generate_plots(wildcards):
    analysis_id=ANALYSES_IDS[0]
    if config["plot"]["generate_plots"]:
        return f"results/{analysis_id}/7_ANI/.gggenomes_done"
    else:
        return f"results/{analysis_id}/prefinal_table.txt"

def gather_gggenomes(wildcards):
    checkpoint_output = checkpoints.prepare_gggenomes_data.get(**wildcards).output[0]
    gggenomes = expand("results/{analysis_id}/7_ANI/2_plot/{gggdata}.png",
                    analysis_id=ANALYSES_IDS[0],
                    gggdata=glob_wildcards(f"{checkpoint_output}/{{gggdata}}.blast").gggdata,
                    )

    return gggenomes
