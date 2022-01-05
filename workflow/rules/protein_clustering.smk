rule proteins_clustering:
    input:
        get_prots_files
    output:
        tsv = f"results/{ANALYSIS_ID}" + "/6_clustering/table_clustering.tsv",
        tmp_dir = temp(directory(f"results/{ANALYSIS_ID}" + "/6_clustering/tmp"))
    conda:
        "../envs/clustering.yaml"
    params:
        ref_faa   = "resources/crassus_dependencies/all_reference_proteins.faa", # TO SET
        prots_faa = f"results/{ANALYSIS_ID}" + "/6_clustering/db/all_proteins.faa",
        prots_db_dir = f"results/{ANALYSIS_ID}" + "/6_clustering/db",
        prots_db   = f"results/{ANALYSIS_ID}" + "/6_clustering/db/all_proteins",
        out_prefix = f"results/{ANALYSIS_ID}" + "/6_clustering/clustering"
    threads: 6
    log:
        db    = f"logs/{ANALYSIS_ID}" + "/clustering/db.log",
        clust = f"logs/{ANALYSIS_ID}" + "/clustering/clustering.log",
        tsv   = f"logs/{ANALYSIS_ID}" + "/clustering/to_table.log"
    shell:
        """
        mkdir -p {params.prots_db_dir} ;
        cat {input} {params.ref_faa} > {params.prots_faa} ;
        mmseqs createdb {params.prots_faa} {params.prots_db} >> {log.db};
        mmseqs cluster {params.prots_db} {params.out_prefix} {output.tmp_dir} \
            -c 0.8 --min-seq-id 0 -s 7.5 --cov-mode 1 --threads {threads} \
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
        presabs = f"results/{ANALYSIS_ID}" + "/6_clustering/presabs_matrix.txt",
        shared  = f"results/{ANALYSIS_ID}" + "/6_clustering/shared_content_matrix.txt",
        nprots  = f"results/{ANALYSIS_ID}" + "/6_clustering/nprots_cluster.txt",
        table_clustering_ids = f"results/{ANALYSIS_ID}" + "/6_clustering/table_clustering_ids.tsv"
    params:
        taxonomy = "resources/crassus_dependencies/reference_taxonomy.txt" # TO SET
    threads: 999
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/calculate_shared_content.py"

rule protein_content_taxa:
    input:
        matrix_shared = rules.calculate_shared_prots.output.shared,
        markers_table = "results/5_phylogenies/taxonomic_classification.txt"
    output:
        "results/6_clustering/shared_content_taxonomy.txt"
    params:
        taxonomy = "resources/crass_taxonomy.txt"
    conda:
        "../../envs/phylogenies.yaml"
    script:
        "../../scripts/taxonomy_from_clustering.py"
