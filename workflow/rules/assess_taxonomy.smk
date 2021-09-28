rule run_DTR_blast:
    input:
        "results/3_crass_contigs/{contig}.fasta"
    output:
        temp(multiext("results/3_crass_contigs/dtr_blast/{contig}.dtr_", "fasta", "db.ndb", "db.nhr", "db.nin", "db.not", "db.nsq", "db.ntf", "db.nto")),
        blast = "results/3_crass_contigs/dtr_blast/{contig}.dtr_blast",
        done = temp("results/3_crass_contigs/dtr_blast/{contig}.dtr_blast_done")
        #done = "results/3_crass_contigs/dtr_blast/{contig}.dtr_blast_done"
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
        "results/5_phylogenies/2_trees/taxonomic_classification.txt"
    params:
        taxonomy = "resources/crass_taxonomy.txt"
    conda:
        "../../envs/phylogenies.yaml"
    script:
        "../../scripts/get_taxonomy_from_trees.py"

rule assess_completenes:
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

rule protein_content_taxa:
    input:
        matrix_shared = rules.calculate_shared_prots.output.shared,
        markers_table  = rules.assess_completenes.output
    output:
        "results/6_clustering/shared_content_taxonomy.txt"
    params:
        taxonomy = "resources/crass_taxonomy.txt"
    conda:
        "../../envs/phylogenies.yaml"
    script:
        "../../scripts/taxonomy_from_clustering.py"

# rule aggregate_taxa:
#     input:
#         #markers_trees = gather_trees,
#         markers_trees = [
#                          "results/5_phylogenies/2_trees/TerL_trimmed.nwk",
#                          "results/5_phylogenies/2_trees/MCP_trimmed.nwk",
#                          "results/5_phylogenies/2_trees/portal_trimmed.nwk",
#                         ],
#         taxa_shared = rules.protein_content_taxa.output
#     output:
#         "results/5_phylogenies/taxonomic_classification_completeness_protshared.txt"
#     params:
#         taxonomy = "resources/crass_taxonomy.txt"
#     conda:
#         "../../envs/phylogenies.yaml"
#     script:
#         "../../scripts/protshared_taxa_to_trees.py"
