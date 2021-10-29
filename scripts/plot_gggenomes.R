library(gggenomes)
library(dplyr)


process_contig_alignment <- function(megablast_file, output_file) {
  megablast <- read.csv(megablast_file, sep = '\t', header = F)
  # create link object
  links <- tibble(seq_id = megablast$V1,
                  start = megablast$V8,
                  end = megablast$V9,
                  seq_id2 = megablast$V2,
                  start2 = megablast$V11,
                  end2 = megablast$V12
  )
  
  # get query and ref contigs
  query <- megablast$V1[1]
  ref   <- megablast$V2[1]
  # get their lengths from blast file too
  query_len <- megablast$V7[1]
  ref_len <- megablast$V10[1]
  # create the genomes track
  genomes <- tibble(seq_id = c(query, ref),
                    length = c(query_len, ref_len)
  )
  genomes
  
  # get query table file
  #qtable_path <- paste(snakemake@params[['contigs_tables_dir']], query, '_prod-', sep='') 
  qtable_path <- paste('results/4_ORF/2_functional_annot_tables/', query, '_tbl-*', sep='') 
  qtable_file <- Sys.glob(qtable_path)[1]
  print(qtable_file)
  
  # get ref table file. Start by checking if it is one of the reference genomes
  rtable_path <- paste('resources/genomes/', ref, '.table', sep='')
  rtable_file <- Sys.glob(rtable_path)[1]
  if (is.na(rtable_file)) {
    rtable_path <- paste('results/4_ORF/2_functional_annot_tables/', ref, '_tbl-*', sep='')
    rtable_file <- Sys.glob(rtable_path)[1]
  }
  print(rtable_file)
  # read genome tables 
  qtable <- read.csv(qtable_file, sep = '\t', header=T)
  rtable <- read.csv(rtable_file, sep = '\t', header=T)
  tables <- rbind(qtable, rtable)
  # create the genes object
  genes <- tibble(seq_id = tables$genome,
                  start = tables$start,
                  end = tables$end,
  )
  
  # check that query and ref are actually different. Otherwise links track will crash
  if (query != ref) {
    p <- gggenomes(genes=genes, seqs=genomes, links=links) + 
      geom_seq() +         # draw contig/chromosome lines
      #geom_bin_label() +   # label each sequence 
      geom_seq_label(nudge_y=-.05) +
      geom_gene() +        # draw genes as arrow
      geom_link() 
    
    # save the plot to png
    ggsave(output_file, plot=p, width=8, height=4)
  } else {
    fileConn<-file(output_file)
    writeLines("self", fileConn)
    close(fileConn) 
  }
}



### main ###

megablast_file <- snakemake@input[['megablast']]
print(megablast_file)
output_file <- snakemake@output[[1]]

# check size of the megablast file. Write empty mock file as output if the former is empty
if (file.size(megablast_file) != 0) {
  process_contig_alignment(megablast_file, output_file)
} else {
  file.create(output_file)
}
  
  
  
  
## geom_seq_label(nudge_y=-.05) 
  


