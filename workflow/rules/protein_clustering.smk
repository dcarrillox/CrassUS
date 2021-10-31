rule proteins_clustering:
    input:
        get_prots_files
    output:
        tsv = temp("results/6_clustering/table_clustering.tsv"),
        tmp_dir = temp(directory("results/6_clustering/tmp"))
    #container: "library://dcarrillo/default/crassus:0.1"
    conda:
        "../../envs/clustering.yaml"
    params:
        prots_faa = "results/6_clustering/db/all_proteins.faa",
        prots_db  = "results/6_clustering/db/all_proteins",
        out_prefx = "results/6_clustering/clustering"
    threads: 4
    log:
        db = "logs/clustering/db.log",
        clust = "logs/clustering/clustering.log",
        tsv = "logs/clustering/to_table.log"
    shell:
        """
        mkdir -p results/6_clustering/db ;
        cat {input} resources/all_crass_proteins.faa > {params.prots_faa} ;
        mmseqs createdb {params.prots_faa} {params.prots_db} >> {log.db};
        mmseqs cluster {params.prots_db} {params.out_prefx} {output.tmp_dir} \
            -c 0 --threads {threads} -s 6 --cluster-steps 4 --cluster-reassign >> {log.clust} ;
        mmseqs createtsv {params.prots_db} {params.prots_db} {params.out_prefx} \
            {output.tsv} >> {log.tsv} ;
        rm -rf results/6_clustering/clustering* results/6_clustering/db
        """

rule calculate_shared_prots:
    input:
        prots_files = get_prots_files,
        tsv = rules.proteins_clustering.output.tsv
    output:
        presabs = "results/6_clustering/presabs_matrix.txt",
        shared  = "results/6_clustering/shared_content_matrix.txt",
        nprots  = temp("results/6_clustering/nprots_cluster.txt")
    params:
        taxonomy = "resources/crass_taxonomy.txt"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/calculate_shared_content.py"

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
