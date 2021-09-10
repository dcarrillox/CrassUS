rule get_marker_proteins:
    input:
        hmmtxt = rules.annotate_proteins_best_coding.output.outfile,
        faa = "results/4_prodigal/best_coding/{prots}.faa",
        gff = "results/4_prodigal/best_coding/{prots}.gff"
    output:
        faa = "results/5_phylogenies/markers_scan/{prots}_markers.faa",
        summary = "results/5_phylogenies/markers_scan/{prots}_markers.summary"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/get_marker_genes.py"

checkpoint summarize_markers:
    input:
        get_markers_files
    output:
        summary = "results/5_phylogenies/markers.summary",
        faa_dir = directory("results/5_phylogenies/markers_faa")
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/summarize_markers.py"

rule multiple_sequence_alignment:
    input:
        found = "results/5_phylogenies/markers_faa/{marker}.faa",
        ref   = "resources/MSAs/{marker}_crassphage_reference.mafft-einsi"
    output:
        "results/5_phylogenies/msa/{marker}.msa"
    threads: 4
    conda:
        "../../envs/phylogenies.yaml"
    shell:
        """
        mafft --add {input.found} --thread {threads} {input.ref} > {output}
        """

rule msa_trimming:
    input:
        "results/5_phylogenies/msa/{marker}.msa"
    output:
        "results/5_phylogenies/msa/{marker}_trimmed.msa"
    conda:
        "../../envs/phylogenies.yaml"
    shell:
        "trimal -in {input} -out {output} -gt 0.9"

rule make_trees:
    input:
        "results/5_phylogenies/msa/{marker}_trimmed.msa"
    output:
        "results/5_phylogenies/tree/{marker}_trimmed.nwk"
    conda:
        "../../envs/phylogenies.yaml"
    shell:
        "fasttree {input} > {output}"

rule measure_leaves_distances:
    input:
        "results/5_phylogenies/tree/{marker}_trimmed.nwk"
    output:
        "results/5_phylogenies/tree/{marker}_trimmed.dist"
    params:
        taxonomy = "resources/crass_taxonomy.txt"
    conda:
        "../../envs/phylogenies.yaml"
    script:
        "../../scripts/measure_leaves_distances.py"
