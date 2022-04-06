rule get_marker_proteins:
    input:
        hmmtxt = rules.annotate_proteins_best_coding.output.outfile,
        faa = "results/{analysis_id}/4_ORF/1_best_coding/{prots}.faa",
        gff = "results/{analysis_id}/4_ORF/1_best_coding/{prots}.gff"
    output:
        faa = "results/{analysis_id}/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.faa",
        summary = "results/{analysis_id}/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.summary"
    params:
        markers_ids = "resources/CrassUS_db/marker_profiles/profiles_length.txt"
    log:
        "logs/{analysis_id}/hmmscan/markers/{prots}_get_markers.log"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/get_marker_genes.py"

rule check_markers:
    input:
        faa = rules.get_marker_proteins.output.faa,
        profile = "resources/CrassUS_db/marker_profiles/custom_yutin_markers.hmm"
    output:
        outfile = "results/{analysis_id}/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.hmmtxt",
        domtblout = "results/{analysis_id}/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.domtxt"
    log:
        "logs/{analysis_id}/hmmscan/markers/{prots}_check_markers.log"
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
            hmmscan --cpu {threads} -o {output.outfile} --domtblout {output.domtblout} {input.profile} {input.faa} 2>{log}
        fi
        '''

checkpoint summarize_markers:
    input:
        get_markers_files
    output:
        summary = "results/{analysis_id}/5_phylogenies/markers_summary.txt",
        coverages = "results/{analysis_id}/5_phylogenies/markers_coverages.txt",
        faa_dir = directory("results/{analysis_id}/5_phylogenies/0_marker_genes/1_final")
    params:
        profiles_length = "resources/CrassUS_db/marker_profiles/profiles_length.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/summarize_markers.py"

rule multiple_sequence_alignment:
    input:
        found = "results/{analysis_id}/5_phylogenies/0_marker_genes/1_final/{marker}.faa",
        ref   = "resources/CrassUS_db/MSAs/reference_{marker}.mafft-einsi"
    output:
        "results/{analysis_id}/5_phylogenies/1_MSAs/{marker}.msa"
    params:
        mode = "mafft-einsi" if config["mafft_msa"]["einsi"] else "mafft-fftnsi --maxiterate 1000"
    threads: 10
    conda:
        "../envs/phylogenies.yaml"
    log:
        "logs/{analysis_id}/msa_alignment/{marker}.log"
    shell:
        '''
        nseqs=$(grep -c ">" {input.found})
        if [[ $nseqs -gt 2000 ]]
        then
            echo "More than 2K sequences, running MAFTT AUTO instead..."
            mafft --auto --quiet --add {input.found} --thread {threads} {input.ref} > {output}
        else
            {params.mode} --quiet --add {input.found} --thread {threads} {input.ref} > {output}
        fi
        '''

rule msa_trimming:
    input:
        rules.multiple_sequence_alignment.output
    output:
        "results/{analysis_id}/5_phylogenies/1_MSAs/{marker}_trimmed.msa"
    conda:
        "../envs/phylogenies.yaml"
    log:
        "logs/{analysis_id}/msa_trimming/{marker}.log"
    shell:
        "trimal -in {input} -out {output} -gt 0.9 &> {log}"

rule make_trees:
    input:
        rules.msa_trimming.output
    output:
        "results/{analysis_id}/5_phylogenies/2_trees/{marker}_trimmed.nwk"
    conda:
        "../envs/phylogenies.yaml"
    log:
        "logs/{analysis_id}/trees/{marker}.log"
    shell:
        "fasttree -quiet -nopr -log {log} {input} > {output}"

rule parse_trees:
    input:
        markers_trees = gather_trees,
        markers_summary = "results/{analysis_id}/5_phylogenies/markers_summary.txt"
    output:
        "results/{analysis_id}/5_phylogenies/markers_classification.txt",
    params:
        taxonomy = "resources/CrassUS_db/reference_taxonomy_subfamily.txt"
    conda:
        "../envs/phylogenies.yaml"
    script:
        "../scripts/get_taxonomy_trees.py"
