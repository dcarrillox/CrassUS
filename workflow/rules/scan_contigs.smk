rule transeq:
    input:
        "results/1_assembly/1_scaff20k/{sample}.fasta"
    output:
        "results/1_assembly/2_six_frames/{sample}.faa",
    shell:
        '''
        size=$(stat --printf="%s" {input})
        if [[ $size -eq 0 ]]
        then
            touch {output}
        else
            transeq -frame 6 -table 11 -clean -sequence {input} -outseq {output}
        fi
        '''

rule hmmsearch_six_frames:
    input:
        "results/1_assembly/2_six_frames/{sample}.faa"
    output:
        txt = "results/2_profiles_scan/{sample}.hmmtxt",
        dom = "results/2_profiles_scan/{sample}.domtxt"
    container: "library://dcarrillo/default/crassus:0.1"
    threads: 2
    shell:
        '''
        size=$(stat --printf="%s" {input})
        if [[ $size -eq 0 ]]
        then
            touch {output}
        else
            hmmsearch --cpu {threads} -o {output.txt} --domtblout {output.dom} --noali --notextw --acc /data/profiles/crass_conserved_genes.hmm {input}
        fi
        '''

rule parse_hmmsearch_six_frames:
    input:
        expand("results/2_profiles_scan/{sample}.hmmtxt", sample=sample_sheet.index)
    output:
        hits_table = "results/2_profiles_scan/profiles_hits.txt",
        contigs    = "results/2_profiles_scan/matching_contigs.txt"
    script:
        "../../scripts/parse_hmmsearch.py"

checkpoint get_matching_contigs:
    input:
        "results/2_profiles_scan/matching_contigs.txt"
    output:
        directory("results/3_contigs/0_contigs")
    params:
        assemblies_dir = "results/1_assembly/1_scaff20k",
        contigs_dir =    "results/3_contigs/0_contigs/"
    script:
        "../../scripts/get_matching_contigs.py"
