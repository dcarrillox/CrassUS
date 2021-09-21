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
    threads: 1
    params:
        length = lambda wildcards, input: input[0].replace(".fasta", "").split("_")[-1]
    log:
        log = "results/4_prodigal/all_codings/{contig}_prod-11.log",
        tga_log = "results/4_prodigal/all_codings/{contig}_prod-TGA.log",
        tag_log = "results/4_prodigal/all_codings/{contig}_prod-TAG.log",
    shell:
        '''
        length={params.length}
        if [[ $length -ge 20000 ]]
        then
            mode="-g 11"
        else
            mode="-p meta"
        fi

        resources/prodigal ${{mode}} -f gff -a {output.faa} -o {output.gff} -i {input} &> {log.log} ;
        resources/prodigal ${{mode}} -TGA W -f gff -a {output.tga_faa} -o {output.tga_gff} -i {input} &> {log.tga_log} ;
        resources/prodigal ${{mode}} -TAG Q -f gff -a {output.tag_faa} -o {output.tag_gff} -i {input} &> {log.tag_log}
        '''

rule coding_density:
    input:
        aggregate_densities
    output:
        "results/4_prodigal/coding_summary.txt"
    params:
        faa_dir = "results/4_prodigal/all_codings"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/summarize_coding_density.py"

checkpoint pick_best_coding:
    input:
        "results/4_prodigal/coding_summary.txt"
    output:
        directory("results/4_prodigal/best_coding")
    params:
        raw_dir = "results/4_prodigal/all_codings"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/pick_best_coding.py"

rule annotate_proteins_best_coding:
    input:
        fasta = "results/4_prodigal/best_coding/{prots}.faa",
        profile = "resources/yutin_2021/all_profiles/all_yutin_profiles.hmm.h3f"
    output:
        outfile = "results/4_prodigal/best_coding/functional_annot/{prots}.hmmtxt",
        domtblout = "results/4_prodigal/best_coding/functional_annot/{prots}.domtxt"
    log:
        "logs/yutin_annot/{prots}.log"
    params:
        evalue_threshold=0.001,
        # if bitscore threshold provided, hmmsearch will use that instead
        #score_threshold=50,
        extra="",
    threads: 2
    wrapper:
        "0.78.0/bio/hmmer/hmmsearch"

rule create_genome_table:
    input:
        gff = "results/4_prodigal/best_coding/{prots}.gff",
        faa = "results/4_prodigal/best_coding/{prots}.faa",
        hmmtxt = "results/4_prodigal/best_coding/functional_annot/{prots}.hmmtxt"
    output:
        "results/4_prodigal/best_coding/genome_tables/{prots}.table"
    params:
        yutin_names = "resources/yutin_2021/all_profiles/yutin_nicknames_profiles.txt"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/create_genome_table.py"


rule genome_tables_finished:
    input:
        get_genome_tables_finished
    output:
        temp("results/4_prodigal/best_coding/genome_tables/.finished")
    shell:
        "touch {output}"
