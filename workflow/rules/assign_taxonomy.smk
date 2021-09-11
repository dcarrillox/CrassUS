rule run_DTR_blast:
    input:
        "results/3_contigs/0_contigs/{contig}.fasta"
    output:
        temp(multiext("results/3_contigs/0_contigs/{contig}.dtr_", "fasta", "db.ndb", "db.nhr", "db.nin", "db.not", "db.nsq", "db.ntf", "db.nto")),
        blast = "results/3_contigs/0_contigs/{contig}.dtr_blast",
        #done = temp("results/3_contigs/0_contigs/{contig}.dtr_blast_done")
        done = "results/3_contigs/0_contigs/{contig}.dtr_blast_done"
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
        # crassus_contigs = expand("results/3_contigs/0_contigs/{contig}.fasta",
        #                         contig=glob_wildcards("results/3_contigs/0_contigs/{contig}.fasta").contig,
        #                         )
        # crassus_contigs = gather_contigs
        # markers_summary = checkpoints.summarize_markers.output.summary
        markers_summary = "results/5_phylogenies/markers.summary"
    output:
        temp("results/5_phylogenies/tree/taxonomic_classification.txt")
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
        "results/5_phylogenies/tree/taxonomic_classification_completeness.txt"
    params:
        lengths = "resources/taxas_average_length.txt"
    conda:
        "../../envs/phylogenies.yaml"
    script:
        "../../scripts/assess_completeness.py"
