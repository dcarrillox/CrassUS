rule blast_all:
    input:
        gather_genomes_blastall,
    output:
        #temp(multiext("results/7_ANI/0_species/all_genomes_ref_crassus.", "fasta", "db.ndb", "db.nhr", "db.nin", "db.not", "db.nsq", "db.ntf", "db.nto"))
        multiext("results/7_ANI/0_species/all_genomes_ref_crassus.", "fasta", "ndb", "nhr", "nin", "not", "nsq", "ntf", "nto"),
        tsv = "results/7_ANI/0_species/all_genomes_ref_crassus_blast.tsv"
    params:
        fasta_all = "results/7_ANI/0_species/all_genomes_ref_crassus.fasta",
        db = "results/7_ANI/0_species/all_genomes_ref_crassus"
    conda:
        "../../envs/compare_genomes.yaml"
    threads: 4
    shell:
        """
        cat {input} > {params.fasta_all} ;
        makeblastdb -in {params.fasta_all} -dbtype nucl -out {params.db} ;
        blastn -query {params.fasta_all} -db {params.db} -dust no \
        -outfmt '6 std qlen slen' -max_target_seqs 10000 \
        -out {output.tsv} -num_threads {threads}
        """

rule anicalc_species:
    input:
        rules.blast_all.output.tsv
    output:
        "results/7_ANI/0_species/anicalc_results.tsv"
    params:
        script = "scripts/anicalc.py"
    conda:
        "../../envs/compare_genomes.yaml"
    shell:
        "python {params.script} -i {input} -o {output}"


rule aniclust_species:
    input:
        ani = rules.anicalc_species.output,
        fasta_all = "results/7_ANI/0_species/all_genomes_ref_crassus.fasta"
    output:
        "results/7_ANI/0_species/aniclust_species.tsv"
    params:
        script = "scripts/aniclust.py"
    conda:
        "../../envs/compare_genomes.yaml"
    shell:
        "python {params.script} --out {output} --fna {input.fasta_all} --ani {input.ani} --min_ani 95 --min_tcov 85 --min_qcov 0"


rule blast_representatives:
    input:
        aniclust = rules.aniclust_species.output,
        fasta_all = "results/7_ANI/0_species/all_genomes_ref_crassus.fasta"
    output:
        #temp(multiext("results/7_ANI/1_genus/repr_genomes.", "fasta", "db.ndb", "db.nhr", "db.nin", "db.not", "db.nsq", "db.ntf", "db.nto"))
        multiext("results/7_ANI/1_genus/repr_genomes.", "fasta", "ndb", "nhr", "nin", "not", "nsq", "ntf", "nto"),
        tsv = "results/7_ANI/1_genus/repr_genomes_blast.tsv",
    params:
        repr_ids = "results/7_ANI/1_genus/repr_genomes.ids",
        fasta_repr = "results/7_ANI/1_genus/repr_genomes.fasta",
        db = "results/7_ANI/1_genus/repr_genomes"
    conda:
        "../../envs/compare_genomes.yaml"
    threads: 4
    shell:
        """
        cut -f1 {input.aniclust} > {params.repr_ids} ;
        seqtk subseq {input.fasta_all} {params.repr_ids} > {params.fasta_repr} ;
        makeblastdb -in {params.fasta_repr} -dbtype nucl -out {params.db} ;
        blastn -query {params.fasta_repr} -db {params.db} \
        -outfmt '6 std qlen slen' -max_target_seqs 10000 \
        -out {output.tsv} -num_threads {threads}
        """


rule anicalc_genus:
    input:
        rules.blast_representatives.output.tsv
    output:
        "results/7_ANI/1_genus/anicalc_results.tsv"
    params:
        script = "scripts/anicalc.py"
    conda:
        "../../envs/compare_genomes.yaml"
    shell:
        "python {params.script} -i {input} -o {output}"


rule aniclust_genus:
    input:
        ani = rules.anicalc_genus.output,
        fasta_repr = "results/7_ANI/1_genus/repr_genomes.fasta"
    output:
        "results/7_ANI/1_genus/aniclust_genus.tsv"
    params:
        script = "scripts/aniclust.py"
    conda:
        "../../envs/compare_genomes.yaml"
    shell:
        "python {params.script} --out {output} --fna {input.fasta_repr} --ani {input.ani} --min_ani 0 --min_tcov 50 --min_qcov 0"


rule get_sp_gen_clusters:
    input:
        sp =  rules.aniclust_species.output,
        gen = rules.aniclust_genus.output,
        fasta_all = "results/7_ANI/0_species/all_genomes_ref_crassus.fasta"
    output:
        "results/7_ANI/ani_genus_species.txt"
    conda:
        "../../envs/compare_genomes.yaml"
    script:
        "../../scripts/create_species_genus_tables_aniclust.py"


rule ani_assing_taxonomy:
    input:
        assignments = rules.get_sp_gen_clusters.output,
        anicalc = rules.anicalc_species.output
    output:
        "results/7_ANI/ani_genus_species_taxonomy.txt"
    params:
        taxonomy = "resources/crass_taxonomy.txt"
    conda:
        "../../envs/compare_genomes.yaml"
    script:
        "../../scripts/get_taxonomy_ani.py"





# aniclust genera
# "python {params.script} --out {output.genus} --fna {params.fasta_all} --ani {input} --min_ani 0 --min_tcov 50 --min_qcov 0 ; "

# rule parse_aniclust_species:
#     input:
#         rules.aniclust_species.output
#     output:



# rule parse_aniclust:
#     input:
#         genus = rules.aniclust.output.genus,
#         species = rules.aniclust.output.species
#     output:
#         genus = "results/7_ANI/genera_clusters.tsv",
#         species = "results/7_ANI/species_clusters.tsv"
#     params:
#         taxonomy = "resources/crass_taxonomy.txt"
#     script:
#         "../../scripts/parse_aniclust.py"
#
#
# #############################
# rule fastani:
#     input:
#         "results/3_crass_contigs/{contig}.fasta",
#     output:
#         "results/7_ANI/0_all/{contig}.fastani"
#     params:
#         fastani_list = "resources/fastANI_genomes.list"
#     threads: 2
#     conda:
#         "../../envs/compare_genomes.yaml"
#     log:
#         "logs/fastani/{contig}.log"
#     shell:
#         "fastANI -q {input} --rl {params.fastani_list} -k 13 --fragLen 1000 "
#         "--minFraction 0.1 -t {threads} -o {output} &> {log}"
#
# rule pyani:
#     input:
#         ani = rules.fastani.output
#     output:
#         "results/7_ANI/1_most_similar/{contig}/ANIb_alignment_coverage.tab"
#     params:
#         outdir = "results/7_ANI/1_most_similar/{contig}",
#         tmp = directory("results/7_ANI/1_most_similar/{contig}_tmp")
#     conda:
#         "../../envs/compare_genomes.yaml"
#     threads: 4
#     log:
#         "logs/pyani/{contig}.log"
#     script:
#         "../../scripts/run_pyani.py"
#
# rule megablast_genomes:
#     input:
#         rules.pyani.output
#     output:
#         megablast = "results/7_ANI/2_plot/{contig}.megablast",
#         tmp  = temp(directory("results/7_ANI/2_plot/{contig}"))
#     params:
#         contigs_dir = "results/3_crass_contigs",
#         refgenomes_dir = "resources/genomes"
#     conda:
#         "../../envs/compare_genomes.yaml"
#     script:
#         "../../scripts/run_megablast.py"
#
# # rule install_gggenomes:
# #     output:
# #         "results/7_ANI/2_plot/.gggenomes_done"
# #     conda:
# #         "../../envs/plot_genomes.yaml"
# #     script:
# #         "../../scripts/install_gggenomes.R"
# #
# #
# rule plot_gggenomes:
#     input:
#         "results/4_ORF/2_functional_annot_tables/.finished",
#         megablast = rules.megablast_genomes.output.megablast,
#     output:
#         "results/7_ANI/2_plot/{contig}.png"
#     params:
#         contigs_tables_dir = "results/4_ORF/2_functional_annot_tables/",
#         reference_tables_dir = "resources/genomes/"
#     # conda:
#     #     "../../envs/plot_genomes.yaml"
#     # script:
#     #     "../../scripts/plot_gggenomes.R"
#     shell:
#         "touch {output}"
