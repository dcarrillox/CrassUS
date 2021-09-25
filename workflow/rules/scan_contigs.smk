rule translate_six_frames:
    input:
        rules.filter_rename_assemblies.output.fasta
    output:
        "results/2_six_frames/0_faa/{sample}.faa"
    conda:
        "../../envs/utils.yaml"
    log:
        "logs/sixframes_transeq/{sample}.log"
    shell:
        '''
        size=$(stat --printf="%s" {input})
        if [[ $size -eq 0 ]]
        then
            touch {output}
        else
            transeq -frame 6 -table 11 -clean -sequence {input} -outseq {output} &> {log}
        fi
        '''

rule hmmsearch_six_frames:
    input:
        fasta = rules.translate_six_frames.output,
        profile = "resources/yutin_2021/scan_contigs/crass_conserved_genes.hmm.h3f"
    output:
        outfile = "results/2_six_frames/1_screening/{sample}.hmmtxt",
        domtblout = "results/2_six_frames/1_screening/{sample}.domtxt"
    log:
        "logs/hmmsearch/sixframes/{sample}.log"
    params:
        evalue_threshold=0.001,
        # if bitscore threshold provided, hmmsearch will use that instead
        #score_threshold=50,
        extra="",
    threads: 2
    wrapper:
        "0.78.0/bio/hmmer/hmmsearch"

rule parse_hmmsearch_six_frames:
    input:
        expand("results/2_six_frames/1_screening/{sample}.hmmtxt", sample=sample_sheet.index)
    output:
        hits_table = "results/2_six_frames/profiles_hits.txt",
        contigs    = "results/2_six_frames/matching_contigs.txt"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/parse_hmmsearch_six_frames.py"

checkpoint get_matching_contigs:
    input:
        rules.parse_hmmsearch_six_frames.output.contigs
    output:
        directory("results/3_crass_contigs/"),
        fastani_list = "resources/fastANI_genomes.list"
    params:
        assemblies_dir = "results/1_rename",
        contigs_dir =    "results/3_crass_contigs",
        genomes_list = "resources/genomes/genomes.list",
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/get_matching_contigs.py"
