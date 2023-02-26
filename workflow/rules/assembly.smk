rule filter_rename_assemblies:
    input:
        get_raw_assemblies,
        "resources/CrassUS_db/.unzip_done"
    output:
        fasta = "results/{analysis_id}/1_rename/{sample}.fasta",
        table = "results/{analysis_id}/1_rename/{sample}.table"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/filter_rename_assemblies.py"


rule run_DTR_blast:
    input:
        "results/{analysis_id}/3_crass_contigs/{contig}.fasta"
    output:
        temp(multiext("results/{analysis_id}/3_crass_contigs/dtr_blast/{contig}.dtr_", "fasta", "db.ndb", "db.nhr", "db.nin", "db.not", "db.nsq", "db.ntf", "db.nto")),
        blast = "results/{analysis_id}/3_crass_contigs/dtr_blast/{contig}.dtr_blast",
        done = "results/{analysis_id}/3_crass_contigs/dtr_blast/{contig}.dtr_blast_done"
    params:
        tmp_fasta = "results/{analysis_id}/3_crass_contigs/dtr_blast/{contig}.dtr_fasta",
        tmp_db = "results/{analysis_id}/3_crass_contigs/dtr_blast/{contig}.dtr_db"
    conda:
        "../envs/utils.yaml"
    log:
        makedb = "logs/{analysis_id}/dtr/{contig}_makedb.log",
        blast  = "logs/{analysis_id}/dtr/{contig}_blast.log"
    script:
        "../scripts/run_DTR_blast.py"
