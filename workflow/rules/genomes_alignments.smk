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
        "results/7_ANI/0_all/{contig}.log"
    shell:
        "fastANI -q {input} --rl {params.fastani_list} -k 13 --fragLen 1000 "
        "--minFraction 0.1 -t {threads} -o {output} &> {log}"


rule pyani:
    input:
        #fasta = "results/3_crass_contigs/{contig}.fasta",
        #ani   = "results/7_ANI/0_all/{contig}.fastani"
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
        "results/7_ANI/1_most_similar/{contig}/{contig}.log"
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
    script:
        "../../scripts/run_megablast.py"

rule genoplot_genomes:
    input:
        "results/4_ORF/2_functional_annot_tables/.finished",
        megablast = rules.megablast_genomes.output.megablast
    output:
        "results/7_ANI/2_plot/{contig}.png"
    params:
        contigs_tables_dir = "results/4_ORF/2_functional_annot_tables"
    # script:
    #     "../../scripts/genoplotr.R"
    shell:
        "touch {output}"
