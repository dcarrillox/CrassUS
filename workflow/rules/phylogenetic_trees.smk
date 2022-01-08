rule get_marker_proteins:
    input:
        hmmtxt = rules.annotate_proteins_best_coding.output.outfile,
        faa = f"results/{ANALYSIS_ID}" + "/4_ORF/1_best_coding/{prots}.faa",
        gff = f"results/{ANALYSIS_ID}" + "/4_ORF/1_best_coding/{prots}.gff"
    output:
        faa = f"results/{ANALYSIS_ID}" + "/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.faa",
        summary = f"results/{ANALYSIS_ID}" + "/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.summary"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/get_marker_genes.py"

rule check_markers:
    input:
        faa = rules.get_marker_proteins.output.faa,
        profile = "resources/crassus_dependencies/marker_profiles/custom_yutin_markers.hmm"
    output:
        outfile = f"results/{ANALYSIS_ID}" + "/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.hmmtxt",
        domtblout = f"results/{ANALYSIS_ID}" + "/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.domtxt"
    # log:
    #     "logs/hmmscan/markers/{prots}.log"
    conda:
        "../envs/utils.yaml"
    threads: 2
    shell:
        '''
        size=$(stat --printf="%s" {input.faa})
        if [[ $size -eq 0 ]]
        then
            touch {output.outfile} ; touch {output.domtblout}
        else
            hmmscan --cpu {threads} -o {output.outfile} --domtblout {output.domtblout} {input.profile} {input.faa}
        fi
        '''

checkpoint summarize_markers:
    input:
        get_markers_files
    output:
        summary = f"results/{ANALYSIS_ID}" + "/5_phylogenies/markers.summary",
        coverages = f"results/{ANALYSIS_ID}" + "/5_phylogenies/markers.coverages",
        faa_dir = directory(f"results/{ANALYSIS_ID}" + "/5_phylogenies/0_marker_genes/1_final")
    params:
        profiles_length = "resources/crassus_dependencies/marker_profiles/profiles_length.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/summarize_markers.py"

if config["mafft_msa"]["einsi"]:
    rule multiple_sequence_alignment_einsi:
        input:
            found = f"results/{ANALYSIS_ID}" + "/5_phylogenies/0_marker_genes/1_final/{marker}.faa",
            ref   = "resources/crassus_dependencies/MSAs/reference_{marker}.mafft-einsi"
        output:
            f"results/{ANALYSIS_ID}" + "/5_phylogenies/1_MSAs/{marker}.msa"
        threads: 10
        conda:
            "../envs/phylogenies.yaml"
        log:
            f"logs/{ANALYSIS_ID}" + "/msa_alignment/{marker}.log"
        shell:
            '''
            nseqs=$(grep -c ">" {input.found})
            if [[ $nseqs -gt 2000 ]]
            then
                echo "More than 2K sequences, running MAFTT AUTO instead..."
                time mafft --auto --quiet --add {input.found} --thread {threads} {input.ref} > {output}
            else
                time mafft-einsi --quiet --add {input.found} --thread {threads} {input.ref} > {output}
            fi
            '''
            #"time mafft-einsi --quiet --add {input.found} --thread {threads} {input.ref} > {output}"
else:
    rule multiple_sequence_alignment_fftnsi:
        input:
            found = f"results/{ANALYSIS_ID}" + "/5_phylogenies/0_marker_genes/1_final/{marker}.faa",
            ref   = "resources/crassus_dependencies/MSAs/reference_{marker}.mafft-einsi"
        output:
            f"results/{ANALYSIS_ID}" + "/5_phylogenies/1_MSAs/{marker}.msa"
        threads: 10
        conda:
            "../envs/phylogenies.yaml"
        log:
            f"logs/{ANALYSIS_ID}" + "/msa_alignment/{marker}.log"
        shell:
            '''
            nseqs=$(grep -c ">" {input.found})
            if [[ $nseqs -gt 2000 ]]
            then
                echo "More than 2K sequences, running MAFTT AUTO instead..."
                time mafft --auto --quiet --add {input.found} --thread {threads} {input.ref} > {output}
            else
                time mafft-fftnsi --maxiterate 1000 --quiet --add {input.found} --thread {threads} {input.ref} > {output}
            fi
            '''
            #"time mafft-fftnsi --maxiterate 1000 --quiet --add {input.found} --thread {threads} {input.ref} > {output}"

rule msa_trimming:
    input:
        f"results/{ANALYSIS_ID}" + "/5_phylogenies/1_MSAs/{marker}.msa"
    output:
        f"results/{ANALYSIS_ID}" + "/5_phylogenies/1_MSAs/{marker}_trimmed.msa"
    conda:
        "../envs/phylogenies.yaml"
    log:
        f"logs/{ANALYSIS_ID}" + "/msa_trimming/{marker}.log"
    shell:
        "trimal -in {input} -out {output} -gt 0.9 &> {log}"

rule make_trees:
    input:
        rules.msa_trimming.output
    output:
        f"results/{ANALYSIS_ID}" + "/5_phylogenies/2_trees/{marker}_trimmed.nwk"
    conda:
        "../envs/phylogenies.yaml"
    log:
        f"logs/{ANALYSIS_ID}" + "/trees/{marker}.log"
    shell:
        "fasttree -log {log} {input} > {output}"



# rule measure_leaves_distances:
#     input:
#         rules.make_trees.output
#     output:
#         "results/5_phylogenies/2_trees/{marker}_trimmed.dist"
#     params:
#         taxonomy = "resources/crass_taxonomy.txt"
#     conda:
#         "../../envs/phylogenies.yaml"
#     script:
#         "../../scripts/measure_leaves_distances.py"
