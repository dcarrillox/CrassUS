import pandas as pd
import glob, os, sys
from snakemake.utils import validate
from snakemake.utils import min_version


configfile: "config/config.yaml"
validate(config, "../schemas/config.schema.yaml")


###### Parse sample sheet ######
sample_sheet = pd.read_table(config["sample_sheet"], comment='#', dtype=str).set_index("sample_id", drop=False)
validate(sample_sheet, "../schemas/samples.schema.yaml")

# replace underscores in the samples_id by hyphen
sample_sheet.index = sample_sheet.index.str.replace('_','-')
# chere sample ids are unique
if sample_sheet.index.has_duplicates:
    print()
    print("samples identifiers are not unique. Exiting...")
    print()
    sys.exit()

# check all the fields are complete by looking at NaN's presence
if sample_sheet.isnull().values.any():
    print()
    print("Not all the fields are complete in the sample sheet. Exiting...")
    print()
    sys.exit()


# check there is only one analysis_id
ANALYSES_IDS = list(set(sample_sheet.analysis_id.tolist()))
if len(ANALYSES_IDS) > 1:
    print()
    print(f"More than one 'analysis_id' in the sample sheet '{os.path.basename(config['sample_sheet'])}':")
    for id in ANALYSES_IDS:
        print(f"\t- {id}")
    print("Please run only one analysis each time. Exiting...")
    print()
    sys.exit()
else:
    SAMPLES=sample_sheet.index.tolist()
#print(sample_sheet.index) # comment when running --dag


###### Wildcard constraints ######
wildcard_constraints:
    analysis_id="|".join(ANALYSES_IDS),
    sample="|".join(sample_sheet.index),
    sample_transeq="|".join(sample_sheet.index),


########## Unzip reference db ############
rule unzip_dependencies:
    input:
        "resources/crassus_dependencies.tar.gz"
    output:
        db_dir = directory("resources/CrassUS_db"),
        mock = "resources/CrassUS_db/.unzip_done"
    params:
        resources_dir = "resources"
    shell:
        "tar zxf {input} --directory {params.resources_dir} ; touch {output.mock}"



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
                        "results/{analysis_id}/4_ORF/3_functional_annot_tables/{best_coding}.table",
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
    ref_genomes = glob.glob("resources/CrassUS_db/reference_genomes/sequences/*.fasta")
    return crassus_fasta + ref_genomes

def gather_trees(wildcards): # used
    # check the config file to know which markers are requested for the trees
    config_markers = [marker for marker in config["phylogenies"]["markers"] if config["phylogenies"]["markers"][marker]]

    # get the markers given by the checkpoint
    checkpoint_output = checkpoints.summarize_markers.get(**wildcards).output.faa_dir
    checkpoint_markers = glob_wildcards(f"{checkpoint_output}/{{marker}}.faa").marker

    # as final, take the intersection between the given by the checpoint and the requested from the config file
    final_markers = [marker for marker in checkpoint_markers if marker in config_markers]

    tree_files = expand("results/{analysis_id}/5_phylogenies/2_trees/{marker}_trimmed.nwk",
                    analysis_id=ANALYSES_IDS[0],
                    marker=final_markers
                    )

    # ask for iToL files as well to trigger its execution.
    tree_itol_files = expand("results/{analysis_id}/5_phylogenies/3_iToL/{marker}_iToL.nwk",
                    analysis_id=ANALYSES_IDS[0],
                    marker=final_markers
                    )
    annot_itol_files = expand("results/{analysis_id}/5_phylogenies/3_iToL/{marker}_iToL.txt",
                    analysis_id=ANALYSES_IDS[0],
                    marker=final_markers
                    )

    return tree_files + tree_itol_files + annot_itol_files

def gather_dtr(wildcards): # used
    checkpoint_output = checkpoints.get_matching_contigs.get(**wildcards).output[0]
    dtr = expand("results/{analysis_id}/3_crass_contigs/dtr_blast/{contig}.dtr_blast_done",
                    analysis_id=ANALYSES_IDS[0],
                    contig=glob_wildcards(f"{checkpoint_output}/{{contig}}.fasta").contig,
                    )
    return dtr

def generate_plots(wildcards):
    analysis_id=ANALYSES_IDS[0]
    if config["genomes_plot"]["generate_plots"]:
        return f"results/{analysis_id}/7_ANI/.gggenomes_done"
    else:
        return f"results/{analysis_id}/crassus_results.tsv"

def gather_gggenomes(wildcards):
    checkpoint_output = checkpoints.prepare_gggenomes_data.get(**wildcards).output[0]
    gggenomes = expand("results/{analysis_id}/7_ANI/2_plot/{gggdata}.png",
                    analysis_id=ANALYSES_IDS[0],
                    gggdata=glob_wildcards(f"{checkpoint_output}/{{gggdata}}.blast").gggdata,
                    )

    return gggenomes
