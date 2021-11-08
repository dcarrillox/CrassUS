library(gggenomes)
library(RColorBrewer)


a <- read.csv("/home/dani/encode_home/projects/crAssUS/writing/Steigviridae_out_tables.txt", sep="\t", header=T, na.strings = "")

genomes_ids <- unique(a$genome)
length_genomes <- unique(a$length)

genomes <- tibble(seq_id = genomes_ids,
                  length = length_genomes 
)

genes <- tibble(seq_id = a$genome,
                start = a$start,
                end = a$end,
                annot = a$yutin,
                shared= a$shared
)

n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))

p <- gggenomes(seqs=genomes, genes=genes) + 
  geom_seq_label() +
  geom_seq() +
  geom_gene(aes(fill=annot, color=shared)) +
  # scale_fill_manual(values=c(MCP="blue",
  #                            portal="red",
  #                            primase="green",
  #                            TerL="orange"), na.value="white") +
  scale_color_manual(values=c(Jelitoviridae="green",
                             Steigviridae="blue"),
                     na.value = "gray50") +
  scale_fill_manual(values=col_vector, na.value="gray90")

pdf("writing/ggtree/Fig_3B.pdf", family = "ArialMT", width = 15.39, height = 8.45, pointsize = 18)
p
invisible(dev.off())

ggsave("writing/ggtree/repr_all", device="svg", dpi=400, width = 2200, height = 800, units = "px")





















----
a <- read.csv("encode_home/no_borrar/coding_slide/2_genome_table/SRR12557716_11_178316_ALL.table", sep="\t", header=F)
names <- c("prot", "genome", "length", "start", "end", "strand", "partial", "yutin")
colnames(a) <- names


genomes_ids <- unique(a$genome)
length_genomes <- unique(a$length)

genomes <- tibble(seq_id = genomes_ids,
                  length = c(rep(length_genomes, 3)) 
                  )

genes <- tibble(seq_id = a$genome,
                start = a$start,
                end = a$end,
                annot = a$yutin
)


p <- gggenomes(seqs=genomes, genes=genes) + 
  geom_seq_label() +
  geom_seq() +
  geom_gene(aes(fill=annot)) +
  scale_fill_manual(values=c(MCP="blue",
                             portal="red",
                             primase="green",
                             TerL="orange"), na.value="grey90")

ggsave("gggenomes_srr.png", plot=p, width=6, height=3)


## two markers plot
library(gggenomes)
library(dplyr)

a <- read.csv("encode_home/no_borrar/two_markers_slide/genome_tables_final.txt", sep="\t", header=F)
names <- c("prot", "genome", "length", "start", "end", "strand", "partial", "yutin")
colnames(a) <- names


#genomes_ids <- unique(a$genome)
genomes_ids <- c("SRR12557730_786_18766", "Yinda_HP26_trimmed_NODE_1_length_94209_cov_109.230165", "SRR12557730_1401_12531")
#length_genomes <- unique(a$length)
length_genomes <- c(18766, 94209, 12531)
  
genomes <- tibble(seq_id = genomes_ids,
                  length = length_genomes) 


genes <- tibble(seq_id = a$genome,
                start = a$start,
                end = a$end,
                annot = a$yutin
)

blast1 <- read.csv("encode_home/no_borrar/two_markers_slide/blast1.txt", sep="\t", header=F)
l0 <- tibble(seq_id=blast1$V1,
             start=blast1$V7,
             end=blast1$V8,
             seq_id2=blast1$V2,
             start2=blast1$V9,
             end2=blast1$V10)

blast2 <- read.csv("encode_home/no_borrar/two_markers_slide/blast2.txt", sep="\t", header=F)
l1 <- tibble(seq_id=blast2$V1,
             start=blast2$V7,
             end=blast2$V8,
             seq_id2=blast2$V2,
             start2=blast2$V9,
             end2=blast2$V10)


l2 <- bind_rows(l0,l1)


p <- gggenomes(seqs=genomes, genes=genes, links=l2) + 
  geom_seq_label() +
  geom_seq() +
  geom_gene(aes(fill=annot)) +
  scale_fill_manual(values=c(MCP="blue",
                             portal="red",
                             primase="green",
                             TerL="orange"), na.value="grey90") +
  geom_link() +
  theme(legend.position = "none")
p
ggsave("gggenomes_two-markers_nolegend.png", plot=p, width=6, height=3)
