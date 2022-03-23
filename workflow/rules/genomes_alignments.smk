rule blast_all:
    input:
        gather_genomes_blastall,
    output:
        temp(multiext("results/{analysis_id}/7_ANI/0_species/all_genomes_ref_crassus.", "fasta", "ndb", "nhr", "nin", "not", "nsq", "ntf", "nto")),
        tsv = "results/{analysis_id}/7_ANI/0_species/all_genomes_ref_crassus_blast.tsv"
    params:
        fasta_all = "results/{analysis_id}/7_ANI/0_species/all_genomes_ref_crassus.fasta",
        db = "results/{analysis_id}/7_ANI/0_species/all_genomes_ref_crassus"
    conda:
        "../envs/compare_genomes.yaml"
    threads: 4
    shell:
        """
        cat {input} > {params.fasta_all} ;
        makeblastdb -in {params.fasta_all} -dbtype nucl -out {params.db} ;
        blastn -query {params.fasta_all} -db {params.db} -dust no \
        -outfmt '6 std qlen slen' -max_target_seqs 10000 \
        -out {output.tsv} -num_threads {threads} -perc_identity 60
        """

rule anicalc_species:
    input:
        rules.blast_all.output.tsv
    output:
        "results/{analysis_id}/7_ANI/0_species/anicalc_results.tsv"
    params:
        script = "workflow/scripts/anicalc.py"
    conda:
        "../envs/compare_genomes.yaml"
    shell:
        "python {params.script} -i {input} -o {output}"

rule aniclust_species:
    input:
        ani = rules.anicalc_species.output,
        fasta_all = "results/{analysis_id}/7_ANI/0_species/all_genomes_ref_crassus.fasta"
    output:
        "results/{analysis_id}/7_ANI/0_species/aniclust_species.tsv"
    params:
        script = "workflow/scripts/aniclust.py"
    conda:
        "../envs/compare_genomes.yaml"
    shell:
        "python {params.script} --out {output} --fna {input.fasta_all} --ani {input.ani} --min_ani 95 --min_tcov 85 --min_qcov 0"

rule blast_representatives:
    input:
        aniclust = rules.aniclust_species.output,
        fasta_all = "results/{analysis_id}/7_ANI/0_species/all_genomes_ref_crassus.fasta"
    output:
        temp(multiext("results/{analysis_id}/7_ANI/1_genus/repr_genomes.", "fasta", "ndb", "nhr", "nin", "not", "nsq", "ntf", "nto")),
        tsv = "results/{analysis_id}/7_ANI/1_genus/repr_genomes_blast.tsv",
    params:
        repr_ids = "results/{analysis_id}/7_ANI/1_genus/repr_genomes.ids",
        fasta_repr = "results/{analysis_id}/7_ANI/1_genus/repr_genomes.fasta",
        db = "results/{analysis_id}/7_ANI/1_genus/repr_genomes"
    conda:
        "../envs/compare_genomes.yaml"
    threads: 4
    shell:
        """
        cut -f1 {input.aniclust} > {params.repr_ids} ;
        seqtk subseq {input.fasta_all} {params.repr_ids} > {params.fasta_repr} ;
        makeblastdb -in {params.fasta_repr} -dbtype nucl -out {params.db} ;
        blastn -query {params.fasta_repr} -db {params.db} \
        -outfmt '6 std qlen slen' -max_target_seqs 10000 \
        -out {output.tsv} -num_threads {threads}
        """

rule anicalc_genus:
    input:
        rules.blast_representatives.output.tsv
    output:
        "results/{analysis_id}/7_ANI/1_genus/anicalc_results.tsv"
    params:
        script = "workflow/scripts/anicalc.py"
    conda:
        "../envs/compare_genomes.yaml"
    shell:
        "python {params.script} -i {input} -o {output}"

rule aniclust_genus:
    input:
        ani = rules.anicalc_genus.output,
        fasta_repr = "results/{analysis_id}/7_ANI/1_genus/repr_genomes.fasta"
    output:
        "results/{analysis_id}/7_ANI/1_genus/aniclust_genus.tsv"
    params:
        script = "workflow/scripts/aniclust.py"
    conda:
        "../envs/compare_genomes.yaml"
    shell:
        "python {params.script} --out {output} --fna {input.fasta_repr} --ani {input.ani} --min_ani 0 --min_tcov 50 --min_qcov 0"

rule get_sp_gen_clusters:
    input:
        sp =  rules.aniclust_species.output,
        gen = rules.aniclust_genus.output,
        fasta_all = "results/{analysis_id}/7_ANI/0_species/all_genomes_ref_crassus.fasta"
    output:
        "results/{analysis_id}/7_ANI/ani_genus_species.txt"
    conda:
        "../envs/compare_genomes.yaml"
    script:
        "../scripts/create_species_genus_tables_aniclust.py"

rule ani_assign_taxonomy:
    input:
        assignments = rules.get_sp_gen_clusters.output,
        anicalc = rules.anicalc_species.output
    output:
        "results/{analysis_id}/7_ANI/ani_genus_species_taxonomy.txt"
    params:
        taxonomy = "resources/crassus_dependencies/reference_taxonomy_subfamily.txt"
    conda:
        "../envs/compare_genomes.yaml"
    script:
        "../scripts/get_taxonomy_ani.py"

rule install_gggenomes:
    output:
        done = "resources/crassus_dependencies/.gggenomes_install_done"
    conda:
        "../envs/plot_genomes.yaml"
    log:
        "logs/install_gggenomes.log"
    shell:
        "Rscript workflow/scripts/install_gggenomes.R &> {log}"

checkpoint prepare_gggenomes_data:
    input:
        "results/{analysis_id}/4_ORF/2_functional_annot_tables/.finished",
        ani = rules.anicalc_species.output,
        blast = rules.blast_all.output.tsv
    output:
        directory("results/{analysis_id}/7_ANI/2_plot")
    params:
        taxonomy = "resources/crassus_dependencies/reference_taxonomy_subfamily.txt",
        genome_tables_dir = "results/{analysis_id}/4_ORF/2_functional_annot_tables",
        ref_genome_tables_dir = "resources/crassus_dependencies/reference_genomes/tables"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/prepare_gggenomes_data.py"

rule plot_gggenomes:
    input:
        rules.install_gggenomes.output.done,
        blast = "results/{analysis_id}/7_ANI/2_plot/{gggdata}.blast",
        annot = "results/{analysis_id}/7_ANI/2_plot/{gggdata}.annot"
    output:
        "results/{analysis_id}/7_ANI/2_plot/{gggdata}.png"
    params:
        contigs_tables_dir = "results/{analysis_id}/4_ORF/2_functional_annot_tables/",
        reference_tables_dir = "resources/crassus_dependencies/reference_genomes/tables"
    conda:
        "../envs/plot_genomes.yaml"
    script:
        "../scripts/plot_gggenomes.R"

rule gather_gggenomes_plots:
    input:
        gather_gggenomes
    output:
        "results/{analysis_id}/7_ANI/.gggenomes_done"
    shell:
        "touch {output}"
