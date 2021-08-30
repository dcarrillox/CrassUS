prodigal_ext = ["prod-11.gff", "prod-11.faa",
                "prod-TGA.gff", "prod-TGA.faa",
                "prod-TAG.gff", "prod-TAG.faa"]
rule prodigal:
    input:
        "results/3_contigs/0_contigs/{contig}.fasta",
    output:
        gff = "results/4_prodigal/all_codings/{contig}_prod-11.gff",
        faa = "results/4_prodigal/all_codings/{contig}_prod-11.faa",
        tga_gff = "results/4_prodigal/all_codings/{contig}_prod-TGA.gff",
        tga_faa = "results/4_prodigal/all_codings/{contig}_prod-TGA.faa",
        tag_gff = "results/4_prodigal/all_codings/{contig}_prod-TAG.gff",
        tag_faa = "results/4_prodigal/all_codings/{contig}_prod-TAG.faa",
    container: "library://dcarrillo/default/crassus:0.1"
    threads: 1
    log:
        log = "results/4_prodigal/all_codings/{contig}_prod-11.log",
        tga_log = "results/4_prodigal/all_codings/{contig}_prod-TGA.log",
        tag_log = "results/4_prodigal/all_codings/{contig}_prod-TAG.log",
    shell:
        "/software/prodigal -g 11 -f gff -a {output.faa} -o {output.gff} -i {input} &> {log.log} ; "
        "/software/prodigal -g 11 -TGA W -f gff -a {output.tga_faa} -o {output.tga_gff} -i {input} &> {log.tga_log} ; "
        "/software/prodigal -g 11 -TAG Q -f gff -a {output.tag_faa} -o {output.tag_gff} -i {input} &> {log.tag_log} "

rule coding_density:
    input:
        # expand("results/4_prodigal/all_codings/{contig}_{prod_ext}",
        #         contig=glob_wildcards("results/3_contigs/0_contigs/{contig}.fasta").contig,
        #         prod_ext=prodigal_ext
        #       )
        aggregate_densities
    output:
        "results/4_prodigal/coding_summary.txt"
    params:
        faa_dir = "results/4_prodigal/all_codings"
    script:
        "../../scripts/summarize_coding_density.py"

checkpoint pick_best_coding:
    input:
        "results/4_prodigal/coding_summary.txt"
    output:
        directory("results/4_prodigal/best_coding")
    params:
        raw_dir = "results/4_prodigal/all_codings"
    script:
        "../../scripts/pick_best_coding.py"
