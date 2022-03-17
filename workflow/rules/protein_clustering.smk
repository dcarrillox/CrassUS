rule proteins_clustering:
    input:
        get_prots_files
    output:
        tsv = "results/{analysis_id}/6_clustering/table_clustering.tsv",
        tmp_dir = temp(directory("results/{analysis_id}/6_clustering/tmp"))
    conda:
        "../envs/clustering.yaml"
    params:
        ref_faa   = "resources/crassus_dependencies/reference_genomes/all_reference_proteins.faa", # TO SET
        prots_faa = "results/{analysis_id}/6_clustering/db/all_proteins.faa",
        prots_db_dir = "results/{analysis_id}/6_clustering/db",
        prots_db   = "results/{analysis_id}/6_clustering/db/all_proteins",
        out_prefix = "results/{analysis_id}/6_clustering/clustering"
    threads: 6
    log:
        db    = "logs/{analysis_id}/clustering/db.log",
        clust = "logs/{analysis_id}/clustering/clustering.log",
        tsv   = "logs/{analysis_id}/clustering/to_table.log"
    shell:
        """
        mkdir -p {params.prots_db_dir} ;
        cat {input} {params.ref_faa} > {params.prots_faa} ;
        mmseqs createdb {params.prots_faa} {params.prots_db} >> {log.db};
        mmseqs cluster {params.prots_db} {params.out_prefix} {output.tmp_dir} \
            -c 0 --min-seq-id 0 -s 7.5 --cluster-steps 4 --threads {threads} \
            --cluster-reassign >> {log.clust} ;
        mmseqs createtsv {params.prots_db} {params.prots_db} {params.out_prefix} \
            {output.tsv} >> {log.tsv} ;
        rm -rf {params.out_prefix}* {params.prots_db_dir}
        """

rule calculate_shared_prots:
    input:
        prots_files = get_prots_files,
        tsv = rules.proteins_clustering.output.tsv
    output:
        presabs = "results/{analysis_id}/6_clustering/presabs_matrix.txt",
        shared  = "results/{analysis_id}/6_clustering/shared_content_matrix.txt",
        nprots  = "results/{analysis_id}/6_clustering/nprots_cluster.txt",
        table_clustering_ids = "results/{analysis_id}/6_clustering/table_clustering_ids.tsv"
    params:
        taxonomy = "resources/crassus_dependencies/reference_taxonomy_subfamily.txt" # TO SET
    threads: 999
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/calculate_shared_content.py"

rule protein_content_taxa:
    input:
        matrix_shared = rules.calculate_shared_prots.output.shared,
        markers_table = "results/{analysis_id}/5_phylogenies/taxonomic_classification.txt"
    output:
        "results/{analysis_id}/6_clustering/shared_content_taxonomy.txt"
    params:
        taxonomy = "resources/crassus_dependencies/reference_taxonomy_subfamily.txt"
    conda:
        "../envs/phylogenies.yaml"
    script:
        "../scripts/taxonomy_from_clustering.py"
