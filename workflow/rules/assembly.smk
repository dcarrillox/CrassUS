rule filter_rename_assemblies:
    input:
        get_raw_assemblies
    output:
        fasta = "results/1_rename/{sample}.fasta",
        table = "results/1_rename/{sample}.table"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/filter_rename_assemblies.py"
