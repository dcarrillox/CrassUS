rule proteins_clustering:
    input:
        get_prots_files
    output:
        tmp = temp(directory("results/6_clustering/tmp")),
        tsv = "results/6_clustering/clustering.tsv"
    container: "library://dcarrillo/default/crassus:0.1"
    params:
        prots_faa = "results/6_clustering/db/all_proteins.faa",
        prots_db  = "results/6_clustering/db/all_proteins",
        out_prefx = "results/6_clustering/clustering"
    threads: 4
    log:
        db = "results/6_clustering/db/db.log",
        clust = "results/6_clustering/clustering.log",
        tsv = "results/6_clustering/tsv.log"
    shell:
        """
        mkdir -p results/6_clustering/db ;
        cat {input} /data/all_crass_proteins.faa > {params.prots_faa} ;
        mmseqs createdb {params.prots_faa} {params.prots_db} >> {log.db};
        mmseqs cluster {params.prots_db} {params.out_prefx} {output.tmp} \
            -c 0 --threads {threads} -s 6 --cluster-steps 4 --cluster-reassign >> {log.clust} ;
        mmseqs createtsv {params.prots_db} {params.prots_db} {params.out_prefx} \
            {output.tsv} >> {log.tsv} ;
        """

rule calculate_shared_prots:
    input:
        prots_files = get_prots_files,
        tsv = rules.proteins_clustering.output.tsv
    output:
        presabs = "results/6_clustering/presabs_matrix.txt",
        shared  = "results/6_clustering/shared_content_matrix.txt"
    params:
        nprots_cluster = "results/6_clustering/nprots_cluster.txt"
    script:
        "../scripts/calculate_shared_content.py"
