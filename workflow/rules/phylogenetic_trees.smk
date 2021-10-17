rule get_marker_proteins:
    input:
        hmmtxt = rules.annotate_proteins_best_coding.output.outfile,
        faa = "results/4_ORF/1_best_coding/{prots}.faa",
        gff = "results/4_ORF/1_best_coding/{prots}.gff"
    output:
        faa = "results/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.faa",
        summary = "results/5_phylogenies/0_marker_genes/0_contigs/{prots}_markers.summary"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/get_marker_genes.py"

checkpoint summarize_markers:
    input:
        get_markers_files
    output:
        summary = "results/5_phylogenies/markers.summary",
        faa_dir = directory("results/5_phylogenies/0_marker_genes/1_final")
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/summarize_markers.py"

if config["mafft_msa"]["einsi"]:
    rule multiple_sequence_alignment_einsi:
        input:
            found = "results/5_phylogenies/0_marker_genes/1_final/{marker}.faa",
            ref   = "resources/MSAs/{marker}_crassphage_reference.mafft-einsi"
        output:
            "results/5_phylogenies/1_MSAs/{marker}.msa"
        threads: 9999
        conda:
            "../../envs/phylogenies.yaml"
        log:
            "logs/msa_alignment/{marker}.log"
        shell:
            "time mafft-einsi --quiet --add {input.found} --thread {threads} {input.ref} > {output}"
else:
    rule multiple_sequence_alignment_fftnsi:
        input:
            found = "results/5_phylogenies/0_marker_genes/1_final/{marker}.faa",
            ref   = "resources/MSAs/{marker}_crassphage_reference.mafft-einsi"
        output:
            "results/5_phylogenies/1_MSAs/{marker}.msa"
        threads: 9999
        conda:
            "../../envs/phylogenies.yaml"
        log:
            "logs/msa_alignment/{marker}.log"
        shell:
            "time mafft-fftnsi --maxiterate 1000 --quiet --add {input.found} --thread {threads} {input.ref} > {output}"

rule msa_trimming:
    input:
        "results/5_phylogenies/1_MSAs/{marker}.msa"
    output:
        "results/5_phylogenies/1_MSAs/{marker}_trimmed.msa"
    conda:
        "../../envs/phylogenies.yaml"
    log:
        "logs/msa_trimming/{marker}.log"
    shell:
        "trimal -in {input} -out {output} -gt 0.9 &> {log}"

rule make_trees:
    input:
        rules.msa_trimming.output
    output:
        "results/5_phylogenies/2_trees/{marker}_trimmed.nwk"
    conda:
        "../../envs/phylogenies.yaml"
    log:
        "logs/trees/{marker}.log"
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
