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
        done = "results/7_pyani/{contig}/.done"
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

rule tblastx_genomes:
    input:
        rules.pyani.output.done
    output:
        done = "results/8_tblastx/{contig}/.done"
    params:
        outdir = "results/8_tblastx/{contig}"
    script:
        "../../scripts/run_tblastx.py"
