options(warn=-1)
suppressPackageStartupMessages(library(gggenomes, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
library(ggnewscale)
library(RColorBrewer)


plot_genomes <- function(blast_file, annot_file, output_file) {
  # print("lol ->", blast_file, annot_file, output_file)
  blast = read.csv(blast_file, sep = '\t', header = T)

  # check how many genomes have to be shown
  tgenomes = unique(blast$tname)
  length(tgenomes)
  if (length(tgenomes) == 1) {
    genomes_ids = c(blast$qname[1], blast$tname[1])
    genomes_len = c(blast$qlen[1], blast$slen[1])
    links <- tibble(seq_id = blast$qname,
                    start = blast$qstart,
                    end = blast$qend,
                    seq_id2 = blast$tname,
                    start2 = blast$tstart,
                    end2 = blast$tend,
                    perc_id = blast$perc_sim
    )
  } else {
    no_ref_df = blast[blast$group == "noref",]
    ref_df    = blast[blast$group == "ref",]
    genomes_ids = c(no_ref_df$tname[1], blast$qname[1], ref_df$tname[1])
    genomes_len = c(no_ref_df$slen[1], blast$qlen[1], ref_df$slen[1])

    l0 <- tibble(seq_id = no_ref_df$qname,
                 start = no_ref_df$qstart,
                 end = no_ref_df$qend,
                 seq_id2 = no_ref_df$tname,
                 start2 = no_ref_df$tstart,
                 end2 = no_ref_df$tend,
                 perc_id = no_ref_df$perc_sim
    )
    links = bind_rows(l0, tibble(seq_id = ref_df$qname,
                                 start = ref_df$qstart,
                                 end = ref_df$qend,
                                 seq_id2 = ref_df$tname,
                                 start2 = ref_df$tstart,
                                 end2 = ref_df$tend,
                                 perc_id = ref_df$perc_sim
    )
    )
  }

  genomes <- tibble(seq_id = genomes_ids,
                    length = genomes_len)

  annot = read.csv(annot_file, sep="\t", header=T)
  genes <- tibble(seq_id = annot$genome,
                  start = annot$start,
                  end = annot$end,
                  annot = annot$yutin
  )

  p <- gggenomes(seqs=genomes, genes=genes, links=links) +
    geom_seq_label() +
    geom_seq() +
    geom_link(aes(fill=perc_id), size=0.02) +
    scale_fill_gradient2(low="red",
                         mid="orange",
                         high="green",
                         limits=c(60, 100),
                         midpoint = 80,
                         name="%id")  +
    new_scale("fill") +
    geom_gene(aes(fill=annot)) +
    scale_fill_manual(values=c(MCP="#1B9E77",
                               portal="#D95F02",
                               primase="#fcf87f",
                               TerL="#E7298A",
                               DNApB="#66A61E",
                               PolA="#E6AB02",
                               Tstab="#A6761D",
                               Ttub="#666666",
                               tail_fib="#a2c4c9",
                               Others="#7570B3"),
                      na.value="white",
                      name="Function")
    #scale_fill_brewer("Genes", palette="Dark2", na.value="gray98")



  ggsave(output_file, plot=p, width=10, height=5)


}



### main ###

blast_file  <- snakemake@input[['blast']]
annot_file  <- snakemake@input[['annot']]
output_file <- snakemake@output[[1]]


# print("here ->", blast_file, annot_file, output_file)

plot_genomes(blast_file, annot_file, output_file)
