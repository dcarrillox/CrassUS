rule filter_rename_assemblies:
    input:
        get_raw_assemblies
    output:
        fasta = "results/{analysis_id}/1_rename/{sample}.fasta",
        table = "results/{analysis_id}/1_rename/{sample}.table"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/filter_rename_assemblies.py"
