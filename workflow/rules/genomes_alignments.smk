rule fastani:
    input:
        "results/3_crass_contigs/{contig}.fasta",
    output:
        "results/7_ANI/0_all/{contig}.fastani"
    params:
        fastani_list = "resources/fastANI_genomes.list"
    threads: 2
    conda:
        "../../envs/compare_genomes.yaml"
    log:
        "logs/fastani/{contig}.log"
    shell:
        "fastANI -q {input} --rl {params.fastani_list} -k 13 --fragLen 1000 "
        "--minFraction 0.1 -t {threads} -o {output} &> {log}"

rule pyani:
    input:
        ani = rules.fastani.output
    output:
        "results/7_ANI/1_most_similar/{contig}/ANIb_alignment_coverage.tab"
    params:
        outdir = "results/7_ANI/1_most_similar/{contig}",
        tmp = directory("results/7_ANI/1_most_similar/{contig}_tmp")
    conda:
        "../../envs/compare_genomes.yaml"
    threads: 4
    log:
        "logs/pyani/{contig}.log"
    script:
        "../../scripts/run_pyani.py"

rule megablast_genomes:
    input:
        rules.pyani.output
    output:
        megablast = "results/7_ANI/2_plot/{contig}.megablast",
        tmp  = temp(directory("results/7_ANI/2_plot/{contig}"))
    params:
        contigs_dir = "results/3_crass_contigs",
        refgenomes_dir = "resources/genomes"
    conda:
        "../../envs/compare_genomes.yaml"
    script:
        "../../scripts/run_megablast.py"

# rule install_gggenomes:
#     output:
#         "results/7_ANI/2_plot/.gggenomes_done"
#     conda:
#         "../../envs/plot_genomes.yaml"
#     script:
#         "../../scripts/install_gggenomes.R"
#
#
rule plot_gggenomes:
    input:
        "results/4_ORF/2_functional_annot_tables/.finished",
        megablast = rules.megablast_genomes.output.megablast,
    output:
        "results/7_ANI/2_plot/{contig}.png"
    params:
        contigs_tables_dir = "results/4_ORF/2_functional_annot_tables/",
        reference_tables_dir = "resources/genomes/"
    # conda:
    #     "../../envs/plot_genomes.yaml"
    # script:
    #     "../../scripts/plot_gggenomes.R"
    shell:
        "touch {output}"
