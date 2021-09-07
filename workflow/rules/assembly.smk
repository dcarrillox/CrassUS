rule filter_rehead_assemblies:
    input:
        get_raw_assemblies
    output:
        fasta = "results/1_assembly/1_scaff20k/{sample}.fasta",
        table = "results/1_assembly/1_scaff20k/{sample}.table"
    conda:
        "../../envs/utils.yaml"
    script:
        "../../scripts/filter_rehead_assemblies.py"
