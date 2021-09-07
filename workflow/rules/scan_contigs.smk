rule translate_six_frames:
    input:
        "results/1_assembly/1_scaff20k/{sample}.fasta"
    output:
        "results/1_assembly/2_six_frames/{sample}.faa"
    conda:
        "../../envs/utils.yaml"
    shell:
        '''
        size=$(stat --printf="%s" {input})
        if [[ $size -eq 0 ]]
        then
            touch {output}
        else
            transeq -frame 6 -table 11 -clean -sequence {input} -outseq {output}
        fi
        '''

rule hmmsearch_six_frames:
    input:
        fasta = "results/1_assembly/2_six_frames/{sample}.faa",
        profile = "resources/yutin_2021/scan_contigs/crass_conserved_genes.hmm.h3f"
    output:
        outfile = "results/2_profiles_scan/{sample}.hmmtxt",
        domtblout = "results/2_profiles_scan/{sample}.domtxt"
    log:
        "logs/hmmsearch/{sample}.log"
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
        expand("results/2_profiles_scan/{sample}.hmmtxt", sample=sample_sheet.index)
    output:
        hits_table = "results/2_profiles_scan/profiles_hits.txt",
        contigs    = "results/2_profiles_scan/matching_contigs.txt"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/parse_hmmsearch_six_frames.py"

checkpoint get_matching_contigs:
    input:
        "results/2_profiles_scan/matching_contigs.txt"
    output:
        directory("results/3_contigs/0_contigs")
    params:
        assemblies_dir = "results/1_assembly/1_scaff20k",
        contigs_dir =    "results/3_contigs/0_contigs/"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/get_matching_contigs.py"
