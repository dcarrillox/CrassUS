rule fastani:
    input:
        "results/3_contigs/0_contigs/{contig}.fasta",
    output:
        "results/6_fastani/{contig}.fastani"
    params:
        fastani_list = "resources/fastANI_genomes.list"
    threads: 2
    conda:
        "../../envs/compare_genomes.yaml"
    log:
        "results/6_fastani/{contig}.log"
    shell:
        "fastANI -q {input} --rl {params.fastani_list} -k 13 --fragLen 1000 "
        "--minFraction 0.1 -t {threads} -o {output} &> {log}"


rule pyani:
    input:
        #fasta = "results/3_contigs/0_contigs/{contig}.fasta",
        #ani   = "results/6_fastani/{contig}.fastani"
        ani = rules.fastani.output
    output:
        "results/7_pyani/{contig}/ANIb_alignment_coverage.tab"
    params:
        outdir = "results/7_pyani/{contig}",
        tmp = directory("results/7_pyani/{contig}_tmp")
    conda:
        "../../envs/compare_genomes.yaml"
    threads: 4
    log:
        "results/7_pyani/{contig}/{contig}.log"
    script:
        "../../scripts/run_pyani.py"

rule megablast_genomes:
    input:
        rules.pyani.output
    output:
        megablast = "results/8_megablast/{contig}.megablast",
        tmp  = temp(directory("results/8_megablast/{contig}"))
    params:
        contigs_dir = "results/3_contigs/0_contigs",
        refgenomes_dir = "resources/genomes"
    script:
        "../../scripts/run_megablast.py"

rule genoplot_genomes:
    input:
        "results/4_prodigal/best_coding/genome_tables/.finished",
        rules.megablast_genomes.output.megablast
    output:
        "results/9_plots/{contig}.txt"
    shell:
        "echo {input} > {output}"
