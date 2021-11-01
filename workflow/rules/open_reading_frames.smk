prodigal_ext = ["tbl-11.gff", "tbl-11.faa",
                "tbl-TGA.gff", "tbl-TGA.faa",
                "tbl-TAG.gff", "tbl-TAG.faa"]
rule predict_ORF:
    input:
        "results/3_crass_contigs/{contig}.fasta",
    output:
        gff = "results/4_ORF/0_all_codings/{contig}_tbl-11.gff",
        faa = "results/4_ORF/0_all_codings/{contig}_tbl-11.faa",
        tga_gff = "results/4_ORF/0_all_codings/{contig}_tbl-TGA.gff",
        tga_faa = "results/4_ORF/0_all_codings/{contig}_tbl-TGA.faa",
        tag_gff = "results/4_ORF/0_all_codings/{contig}_tbl-TAG.gff",
        tag_faa = "results/4_ORF/0_all_codings/{contig}_tbl-TAG.faa",
    threads: 1
    params:
        length = lambda wildcards, input: input[0].replace(".fasta", "").split("_")[-1]
    log:
        log = "logs/orf_prediction/{contig}_tbl-11.log",
        tga_log = "logs/orf_prediction/{contig}_tbl-TGA.log",
        tag_log = "logs/orf_prediction/{contig}_tbl-TAG.log",
    shell:
        '''
        length={params.length}
        if [[ $length -ge 100000 ]]
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
        "results/4_ORF/coding_summary.txt"
    params:
        faa_dir = "results/4_ORF/all_codings"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/summarize_coding_density.py"

checkpoint pick_best_coding:
    input:
        "results/4_ORF/coding_summary.txt"
    output:
        directory("results/4_ORF/1_best_coding")
    params:
        raw_dir = "results/4_ORF/0_all_codings"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/pick_best_coding.py"

rule annotate_proteins_best_coding:
    input:
        fasta = "results/4_ORF/1_best_coding/{prots}.faa",
        #profile = "resources/yutin_2021/all_profiles/all_yutin_profiles.hmm.h3f"
        profile = "resources/yutin_2021/all_profiles/yutin_and_markers.hmm.h3f"
    output:
        outfile = "results/4_ORF/2_functional_annot/{prots}.hmmtxt",
        domtblout = "results/4_ORF/2_functional_annot/{prots}.domtxt"
    log:
        "logs/hmmsearch/orfs/{prots}.log"
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
        gff = "results/4_ORF/1_best_coding/{prots}.gff",
        faa = "results/4_ORF/1_best_coding/{prots}.faa",
        hmmtxt = rules.annotate_proteins_best_coding.output.outfile,
        prots_clusters = "results/6_clustering/table_clustering_ids.tsv"
    output:
        "results/4_ORF/2_functional_annot_tables/{prots}.table"
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
        #temp("results/4_ORF/2_functional_annot_tables/.finished"),
        "results/4_ORF/2_functional_annot_tables/.finished"
    shell:
        "touch {output}"
