rule get_marker_proteins:
    input:
        hmmtxt = rules.annotate_proteins_best_coding.output.outfile,
        faa = "results/4_prodigal/best_coding/{prots}.faa",
        gff = "results/4_prodigal/best_coding/{prots}.gff"
    output:
        faa = "results/5_phylogenies/markers_scan/{prots}_markers.faa",
        summary = "results/5_phylogenies/markers_scan/{prots}_markers.summary"
    #container: "library://dcarrillo/default/crassus:0.1"
    script:
        "../../scripts/get_marker_genes.py"

rule mafft_terl:
    input:
        get_terl_faa_files
    output:
        query_terl = "results/6_terl_tree/query_terl_seqs.faa",
        mafft_out = "results/6_terl_tree/terl_alignment.fasta"
    container: "library://dcarrillo/default/crassus:0.1"
    threads: 4
    shell:
        """
        cat {input} > {output.query_terl}
        size=$(stat --printf="%s" {output.query_terl})
        if [[ $size -eq 0 ]]
        then
            touch {output.mafft_out}
        else
            mafft --add {output.query_terl} --thread 4 /data/crass_reference_TerL.mafft-einsi > {output.mafft_out}
        fi
        """

rule fasttree_terl:
    input:
        rules.mafft_terl.output.mafft_out
    output:
        "results/6_terl_tree/terl_tree.nwk"
    container: "library://dcarrillo/default/crassus:0.1"
    log:
        "results/6_terl_tree/terl_tree.log"
    shell:
        """
        size=$(stat --printf="%s" {input})
        if [[ $size -eq 0 ]]
        then
            touch {output}
        else
            fasttree -log {log} {input} > {output}
        fi
        """
