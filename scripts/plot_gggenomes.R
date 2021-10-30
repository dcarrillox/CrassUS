suppressPackageStartupMessages(library(gggenomes, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))


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
                    end2 = blast$tend
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
                 end2 = no_ref_df$tend
    )
    links = bind_rows(l0, tibble(seq_id = ref_df$qname,
                                 start = ref_df$qstart,
                                 end = ref_df$qend,
                                 seq_id2 = ref_df$tname,
                                 start2 = ref_df$tstart,
                                 end2 = ref_df$tend
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
    geom_link() +
    geom_gene(aes(fill=annot)) +
    # scale_fill_manual(values=c(MCP="blue",
    #                            portal="red",
    #                            primase="green",
    #                            TerL="orange"), na.value="grey90")
    scale_fill_brewer("Genes", palette="Dark2", na.value="gray98")
  
    ggsave(output_file, plot=p, width=8, height=4)
  
  
  # # get query and ref contigs
  # query <- megablast$V1[1]
  # ref   <- megablast$V2[1]
  # # get their lengths from blast file too
  # query_len <- megablast$V7[1]
  # ref_len <- megablast$V10[1]
  # # create the genomes track
  # genomes <- tibble(seq_id = c(query, ref),
  #                   length = c(query_len, ref_len)
  # )
  # genomes
  # 
  # # get query table file
  # #qtable_path <- paste(snakemake@params[['contigs_tables_dir']], query, '_prod-', sep='') 
  # qtable_path <- paste('results/4_ORF/2_functional_annot_tables/', query, '_tbl-*', sep='') 
  # qtable_file <- Sys.glob(qtable_path)[1]
  # print(qtable_file)
  # 
  # # get ref table file. Start by checking if it is one of the reference genomes
  # rtable_path <- paste('resources/genomes/', ref, '.table', sep='')
  # rtable_file <- Sys.glob(rtable_path)[1]
  # if (is.na(rtable_file)) {
  #   rtable_path <- paste('results/4_ORF/2_functional_annot_tables/', ref, '_tbl-*', sep='')
  #   rtable_file <- Sys.glob(rtable_path)[1]
  # }
  # print(rtable_file)
  # # read genome tables 
  # qtable <- read.csv(qtable_file, sep = '\t', header=T)
  # rtable <- read.csv(rtable_file, sep = '\t', header=T)
  # tables <- rbind(qtable, rtable)
  # # create the genes object
  # genes <- tibble(seq_id = tables$genome,
  #                 start = tables$start,
  #                 end = tables$end,
  # )
  # 
  # # check that query and ref are actually different. Otherwise links track will crash
  # if (query != ref) {
  #   p <- gggenomes(genes=genes, seqs=genomes, links=links) + 
  #     geom_seq() +         # draw contig/chromosome lines
  #     #geom_bin_label() +   # label each sequence 
  #     geom_seq_label(nudge_y=-.05) +
  #     geom_gene() +        # draw genes as arrow
  #     geom_link() 
  #   
  #   # save the plot to png
  #   ggsave(output_file, plot=p, width=8, height=4)
  # } else {
  #   fileConn<-file(output_file)
  #   writeLines("self", fileConn)
  #   close(fileConn) 
  # }
}



### main ###

blast_file  <- snakemake@input[['blast']]
annot_file  <- snakemake@input[['annot']]
output_file <- snakemake@output[[1]]


# print("here ->", blast_file, annot_file, output_file)

plot_genomes(blast_file, annot_file, output_file)

  


