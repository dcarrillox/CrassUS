from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

# env: utils.yaml

# read the contig file and get start and end of the sequence
record = SeqIO.read(snakemake.input[0], "fasta")
start = SeqRecord(Seq(record.seq[:1000]), id=f"{record.id}|start", description="")
end   = SeqRecord(Seq(record.seq[-1000:]), id=f"{record.id}|end", description="")

to_write = [start, end]

with open(snakemake.params.tmp_fasta, "w") as fout:
    SeqIO.write(to_write, fout, "fasta")

# build blast database
os.system(f"makeblastdb -in {snakemake.params.tmp_fasta} -out {snakemake.params.tmp_db} -dbtype nucl -logfile {snakemake.log.makedb}")
# run blast
os.system(f"blastn -query {snakemake.params.tmp_fasta} -db {snakemake.params.tmp_db} \
           -out {snakemake.output.blast} -outfmt 6 -logfile {snakemake.log.blast}")

os.system(f"touch {snakemake.output.done}")
