library(genoPlotR)
library(dplyr)

# read megablast file
megablast_file <- snakemake@input[['megablast']]

# create comparison object
blast_coord <- read.csv(megablast_file, sep = '\t', header = F)
blast_coord_start_end <- select(blast_coord, V8, V9, V11, V12)
colnames(blast_coord_start_end) <- c("start1", "end1", "start2", "end2")
comparison <- list(as.comparison(blast_coord_start_end))

# get query and ref contigs
query <- blast_coord$V1[1]
ref   <- blast_coord$V2[1]

# get query table file
#qtable_path <- paste(snakemake@params[['contigs_tables_dir']], query, '_prod-', sep='') 
qtable_path <- paste('/home/danielc/projects/crAssUS/results/4_prodigal/best_coding/genome_tables/', query, '_prod-*', sep='') 
qtable_file <- Sys.glob(qtable_path)[1]
# get ref table file. Start by checking if it is one of the reference genomes
rtable_path <- paste('/home/danielc/projects/crAssUS/resources/genomes/', ref, '.table', sep='')
rtable_file <- Sys.glob(rtable_path)[1]
if (is.na(rtable_file)) {
  rtable_path <- paste('/home/danielc/projects/crAssUS/results/4_prodigal/best_coding/genome_tables/', ref, '.table', sep='')
  rtable_file <- Sys.glob(rtable_path)[1]
}



## ## ## ## ## 
## PLAYING ##
## ## ## ## ## 
megablast_file = '/home/danielc/projects/crAssUS/results/8_megablast/Baboon36B_169_115010.megablast'
table_1 = '/home/danielc/projects/crAssUS/results/4_prodigal/best_coding/genome_tables/Baboon36B_169_115010_prod-11.table'
table_2 = '/home/danielc/projects/crAssUS/results/4_prodigal/best_coding/genome_tables/Baboon22_95_118645_prod-11.table'

dnaseq1 = read.csv(table_1, sep='\t')
dnaseq1

create_dna_seg <- function(table_file) {
  # read genome table file
  df <- read.csv(table_file, sep = '\t')
  # select specific columns
  df2 <- select(df, 'protein_id', 'start', 'stop', 'strand')
  # add direction column & remove the strand one
  df2$direction <- ifelse(df2$strand == "+", 1, -1)
  df2 <- df2 %>% select(-c('strand'))
  # create columns to add
  col <- rep('black', length(rownames(df)))
  lty <- rep(1, length(rownames(df)))
  lwd <- rep(1, length(rownames(df)))
  pch <- rep(1, length(rownames(df)))
  cex <- rep(1, length(rownames(df)))
  gene_type <- rep('arrows', length(rownames(df)))
  fill <- rep('blue', length(rownames(df)))
  # add columns
  df3 <- cbind(df2, col, lty, lwd, pch, cex, gene_type, fill)
  # change colnames to match what genoplotr expects
  colnames(df3) <- c('name', 'start', 'end', 'strand', 'col', 'lty', 'lwd', 'pch', 'cex', 'gene_type', 'fill')
  seg <- dna_seg(df3)
  return(seg)
}


seg1 <- create_dna_seg(table_1)
seg2 <- create_dna_seg(table_2)
class(seg)
colnames(seg)

segs <- list()
segs[[1]] <- seg1
segs[[2]] <- seg2
plot_gene_map(segs)

# read megablast file
blast_coord <- read.csv(megablast_file, sep = '\t', header = F)
blast_coord_start_end <- select(blast_coord, V8, V9, V11, V12)
colnames(blast_coord_start_end) <- c("start1", "end1", "start2", "end2")
comparison <- list(as.comparison(blast_coord_start_end))

plot_gene_map(segs, comparisons = comparison)
