rule aggregate_signals:
    input:
        phylogenies = rules.parse_trees.output,
        shared_prot = "results/{analysis_id}/6_protein_clustering/shared_content_taxonomy.txt",
        ani_cluster = "results/{analysis_id}/7_ANI/ani_genus_species_taxonomy.txt"
    output:
        "results/{analysis_id}/aggregated_signals.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/aggregate_signals.py"

rule final_table:
    input:
        rules.aggregate_signals.output,
        dtr_blast_done = gather_dtr
    output:
        "results/{analysis_id}/crassus_results.tsv"
    params:
        lengths = "resources/CrassUS_db/taxas_average_length.txt",
        ids_dir = "results/{analysis_id}/1_rename",
        taxonomy = 'resources/CrassUS_db/reference_taxonomy_subfamily.txt'
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/final_table.py"
