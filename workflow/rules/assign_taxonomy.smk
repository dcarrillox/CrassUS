rule run_DTR_blast:
    input:
        "results/3_contigs/0_contigs/{contig}.fasta"
    output:
        temp(multiext("results/3_contigs/0_contigs/{contig}.dtr_", "fasta", "db.ndb", "db.nhr", "db.nin", "db.not", "db.nsq", "db.ntf", "db.nto")),
        blast = "results/3_contigs/0_contigs/{contig}.dtr_blast",
        done = "results/3_contigs/0_contigs/.{contig}.dtr_blast_done",
    params:
        tmp_fasta = "results/3_contigs/0_contigs/{contig}.dtr_fasta",
        tmp_db = "results/3_contigs/0_contigs/{contig}.dtr_db"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/run_DTR_blast.py"


rule parse_trees:
    input:
        markers_trees = gather_trees,
        crassus_contigs = expand("results/3_contigs/0_contigs/{contig}.fasta",
                                contig=glob_wildcards("results/3_contigs/0_contigs/{contig}.fasta").contig,
                                )
    output:
        "results/5_phylogenies/tree/taxonomic_classification.txt"
    params:
        taxonomy = "resources/crass_taxonomy.txt"
    conda:
        "../../envs/phylogenies.yaml"
    script:
        "../../scripts/get_taxonomy_from_trees.py"

# rule assign_taxonomy:
#     input:
#
#
#
#         terl_tree = rules.fasttree_terl.output,
#         shared_c  = rules.calculate_shared_prots.output.shared,
#         found_contigs = expand("results/3_contigs/0_contigs/{contig}.fasta",
#                                 contig=glob_wildcards("results/3_contigs/0_contigs/{contig}.fasta").contig,
#                               )
#     output:
#         "results/taxonomic_classification.txt"
#     params:
#         taxonomy = "resources/terL/crass_taxonomy.txt"
#     #container: "library://dcarrillo/default/crassus:0.1"
#     conda: "../../envs/phylogenies.yaml"
#     script: "../scripts/final_taxonomic_classification.py"
