if in_fastq:
    # run assembly
    rule spades_assembly:
        input:
            get_final_reads
        output:
            "results/1_assembly/0_raw/{sample}/contigs.fasta"
        threads:
            2
        params:
            outdir = "results/1_assembly/0_raw/{sample}",
            input_reads = lambda wildcards: format_spades_input(wildcards)
        log:
            "logs/spades/{sample}.log"
        shell:
            "spades.py -t {threads} -o {params.outdir} "
            "{params.input_reads} &>{log}"

rule filter_rehead_assemblies:
    input:
        get_raw_assemblies
    output:
        fasta = "results/1_assembly/1_scaff20k/{sample}.fasta",
        table = "results/1_assembly/1_scaff20k/{sample}.table"
    script:
        "../../scripts/filter_rehead_assemblies.py"
