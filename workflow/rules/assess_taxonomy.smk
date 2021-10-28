rule run_DTR_blast:
    input:
        "results/3_crass_contigs/{contig}.fasta"
    output:
        temp(multiext("results/3_crass_contigs/dtr_blast/{contig}.dtr_", "fasta", "db.ndb", "db.nhr", "db.nin", "db.not", "db.nsq", "db.ntf", "db.nto")),
        blast = "results/3_crass_contigs/dtr_blast/{contig}.dtr_blast",
        #done = temp("results/3_crass_contigs/dtr_blast/{contig}.dtr_blast_done")
        done = "results/3_crass_contigs/dtr_blast/{contig}.dtr_blast_done"
    params:
        tmp_fasta = "results/3_crass_contigs/dtr_blast/{contig}.dtr_fasta",
        tmp_db = "results/3_crass_contigs/dtr_blast/{contig}.dtr_db"
    conda:
        "../../envs/utils.yaml"
    log:
        makedb = "logs/dtr/{contig}_makedb.log",
        blast  = "logs/dtr/{contig}_blast.log"
    script:
        "../../scripts/run_DTR_blast.py"

rule parse_trees:
    input:
        markers_trees = gather_trees,
        markers_summary = "results/5_phylogenies/markers.summary"
    output:
        #temp("results/5_phylogenies/2_trees/taxonomic_classification.txt")
        "results/5_phylogenies/2_trees/taxonomic_classification.txt",
    params:
        taxonomy = "resources/crass_taxonomy.txt"
    conda:
        "../../envs/phylogenies.yaml"
    script:
        "../../scripts/get_taxonomy_from_trees.py"

rule assess_completeness:
    input:
        taxa_markers = rules.parse_trees.output,
        dtr_blast_done = gather_dtr
    output:
        "results/5_phylogenies/taxonomic_classification_completeness.txt"
    params:
        lengths = "resources/taxas_average_length.txt"
    conda:
        "../../envs/phylogenies.yaml"
    script:
        "../../scripts/assess_completeness.py"


rule aggregate_signals:
    input:
        phylogenies = rules.assess_completeness.output,
        shared_prot = "results/6_clustering/shared_content_taxonomy.txt",
        ani_cluster = "results/7_ANI/ani_genus_species_taxonomy.txt"
    output:
        "results/crassus_results.txt"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/final_table.py"

# rule aggregate_taxa_sources:
#     input:
#         markers_trees = gather_trees,
#         taxa_markers  = rules.assess_completenes.output,
#         taxa_shared = rules.protein_content_taxa.output
#     output:
#         "results/5_phylogenies/taxonomic_classification_completeness_protshared.txt"
#     params:
#         taxonomy = "resources/crass_taxonomy.txt",
#         lengths = "resources/taxas_average_length.txt"
#     conda:
#         "../../envs/phylogenies.yaml"
#     script:
#         "../../scripts/protshared_taxa_to_trees.py"
#
# rule assess_unknown_genomes:
#     input:
#         markers_trees = gather_trees,
#         taxa_assessments = rules.aggregate_taxa_sources.output,
#         genome_tables = rules.genome_tables_finished.output,
#         sharing_percentages = rules.calculate_shared_prots.output.shared
#     output:
#         "results/5_phylogenies/unknown_genomes.txt"
#     params:
#         taxonomy = "resources/crass_taxonomy.txt",
#         geno_tables_dir = "results/4_ORF/2_functional_annot_tables"
#     conda:
#         "../../envs/phylogenies.yaml"
#     script:
#         "../../scripts/assess_unknown_genomes.py"
#
#
# rule assess_new_genera:
#     input:
#         markers_trees = gather_trees,
#         taxa_markers  = rules.aggregate_taxa_sources.output,
#         taxa_shared = rules.protein_content_taxa.output
#     output:
#         "results/5_phylogenies/taxonomic_classification_completeness_protshared_newgen.txt"
#     params:
#         taxonomy = "resources/crass_taxonomy.txt"
#     conda:
#         "../../envs/phylogenies.yaml"
#     script:
#         "../../scripts/assess_new_genera.py"
#
# rule assess_species:
#     input:
#         pyani = get_pyani,
#         taxa_table = rules.assess_new_genera.output
#     output:
#         "results/7_ANI/taxonomic_classification_completeness_protshared_newgen_species.txt"
#     params:
#         taxonomy = "resources/crass_taxonomy.txt"
#     conda:
#         "../../envs/utils.yaml"
#     script:
#         "../../scripts/assess_species.py"
