rule translate_six_frames:
    input:
        rules.filter_rename_assemblies.output.fasta
    output:
        "results/{analysis_id}/2_six_frames/0_faa/{sample}.faa"
    conda:
        "../envs/utils.yaml"
    log:
        "logs/{analysis_id}/sixframes_transeq/{sample}.log"
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
        profile = "resources/CrassUS_db/marker_profiles/custom_yutin_markers.hmm.h3f"
    output:
        outfile = "results/{analysis_id}/2_six_frames/1_screening/{sample}.hmmtxt",
        domtblout = "results/{analysis_id}/2_six_frames/1_screening/{sample}.domtxt"
    log:
        "logs/{analysis_id}/hmmsearch/sixframes/{sample}.log"
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
        expand("results/{analysis_id}/2_six_frames/1_screening/{sample}.hmmtxt",
               sample=sample_sheet.index, analysis_id=ANALYSES_IDS)
    output:
        hits_table = "results/{analysis_id}/2_six_frames/profiles_hits.txt",
        contigs    = "results/{analysis_id}/2_six_frames/matching_contigs.txt"
    params:
        profiles_names = "resources/CrassUS_db/functional_annot/nicknames.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/parse_hmmsearch_six_frames.py"

checkpoint get_matching_contigs:
    input:
        rules.parse_hmmsearch_six_frames.output.contigs
    output:
        directory("results/{analysis_id}/3_crass_contigs/"),
    params:
        assemblies_dir = "results/{analysis_id}/1_rename",
        contigs_dir =    "results/{analysis_id}/3_crass_contigs",
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/get_matching_contigs.py"
