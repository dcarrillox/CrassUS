if config["fastq_processing"]:
    rule qc_fastq:
        input: "results/qc/multiqc_trim.html",

#################################

# fastq files as input
if in_fastq:
    # if preprocessing
    if config["fastq_processing"]:

        rule fastqc:
            input:
                unpack(get_fastq),
            output:
                html="results/qc/{sample}.html",
                zip="results/qc/{sample}_fastqc.zip",
            log:
                "logs/fastqc/{sample}.log",
            wrapper:
                "0.74.0/bio/fastqc"

        rule multiqc:
            input:
                expand(["results/qc/{sample}_fastqc.zip"],
                        sample=sample_sheet[sample_sheet['type']=="fastq"]["sample_id"])
            output:
                report(
                    "results/qc/multiqc.html",
                    caption="../report/multiqc.rst",
                    category="Quality control",
                ),
            log:
                "logs/multiqc.log",
            wrapper:
                "0.74.0/bio/multiqc"

        rule trim_reads_se:
            input:
                unpack(get_fastq),
            output:
                "results/trimmed/{sample}.fastq",
            params:
                **config["params"]["trimmomatic"]["se"],
                extra="",
            log:
                "logs/trimmomatic/{sample}.log",
            wrapper:
                "0.74.0/bio/trimmomatic/se"

        rule trim_reads_pe:
            input:
                unpack(get_fastq),
            output:
                r1="results/trimmed/{sample}_R1.fastq",
                r2="results/trimmed/{sample}_R2.fastq",
                r1_unpaired=temp("results/trimmed/{sample}_R1.unpaired.fastq"),
                r2_unpaired=temp("results/trimmed/{sample}_R2.unpaired.fastq"),
                trimlog="results/trimmed/{sample}.trimlog.txt",
            params:
                **config["params"]["trimmomatic"]["pe"],
                extra=lambda w, output: "-trimlog {}".format(output.trimlog),
                #extra=""
            log:
                "logs/trimmomatic/{sample}.log",
            wrapper:
                "0.74.0/bio/trimmomatic/pe"

        rule fastqc_trim:
            input:
                unpack(get_trimmed_reads)
            output:
                html="results/qc/{sample}_trim.html",
                zip="results/qc/{sample}_trim_fastqc.zip",
            log:
                "logs/fastqc/{sample}_trim.log",
            wrapper:
                "0.74.0/bio/fastqc"

        rule multiqc_trim:
            input:
                expand(["results/qc/{sample}_trim_fastqc.zip"],
                        sample=sample_sheet[sample_sheet['type']=="fastq"]["sample_id"])
            output:
                report(
                    "results/qc/multiqc_trim.html",
                    caption="../report/multiqc.rst",
                    category="Quality control",
                ),
            log:
                "logs/multiqc_trim.log",
            wrapper:
                "0.74.0/bio/multiqc"
