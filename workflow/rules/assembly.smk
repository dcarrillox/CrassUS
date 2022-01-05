rule filter_rename_assemblies:
    input:
        get_raw_assemblies
    output:
        fasta = f"results/{ANALYSIS_ID}" + "/1_rename/{sample}.fasta",
        table = f"results/{ANALYSIS_ID}" + "/1_rename/{sample}.table"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/filter_rename_assemblies.py"
