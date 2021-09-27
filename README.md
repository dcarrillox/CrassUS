# crAssUS

Below are detailed the necessary steps to run crAssUS on your set of contigs. 

1) Install Snakemake

You will need `snakemake` to run the pipeline. To ensure the latest version is installed, create a separate environment for it with conda.
```
$ conda create -c bioconda -c conda-forge --name crassus_dev snakemake

$ conda activate crassus_dev
```

2) Clone this repo in your machine

```
$ git clone https://github.com/dcarrillox/crAssUS.git
```

3) Copy the `resources` folder database from my machine.

```
$ cp -r /home/danielc/projects/crAssUS/resources/genomes/ crAssUS/resources/ #yes to overwrite "genomes.list"
```

4) Fill `config/samples.tsv` with the paths to your metagenomes. It is a TAB delimited file, just change the paths and sample_ids:

```
sample_id	fasta
Baboon19	/home/danielc/danielc2/ancient_microbiomes/results/2_assembly/Baboon19/Baboon19.contigs.fa
Chimp20	/home/danielc/danielc2/ancient_microbiomes/results/2_assembly/Chimp20/Chimp20.contigs.fa
Dunia_s_4	/home/danielc/danielc2/ancient_microbiomes/results/2_assembly/Dunia_s_4/Dunia_s_4.contigs.fa
ERR1094777	/home/danielc/danielc2/ancient_microbiomes/results/2_assembly/ERR1094777/ERR1094777.contigs.fa
```

5) Run the pipeline from the `crAssUS` directory. A dry run first

```
$ snakemake -j10 --use-conda -p -n
```

If everything went well, run the pipeline. Adjust the number of CPUs (10 in this case) to your needs:

```
`$ snakemake -j10 --use-conda -p
```

Hopefully it will crAs(s)h. 
