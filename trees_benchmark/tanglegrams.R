library(ape)
library(dendextend)
library(phytools)
library(dplyr)
library(phylogram)

setwd("projects/crAssUS/")
#setwd("~/encode_home/projects/crAssUS/")

tree1 <- read.tree(file = "mafft_test/gold_standard_famcollpsd.nwk")
tree1 <- root(tree1, 'NC_021803|775|90')
tree2 <- read.tree(file = "mafft_test/add/einsi/terl_ref_distant_einsi_trim09_fasttree_famcollpsd.nwk")
tree2 <- root(tree2, 'NC_021803|775|90')
dnd1 <- as.dendrogram.phylo(tree1)
dnd2 <- as.dendrogram.phylo(tree2)

c <-cophylo(tree1, tree2)
plot(c)

dndlist <- dendextend::dendlist(dnd1, dnd2)
dendextend::tanglegram(dndlist, 
                       sort=TRUE,
                       lab.cex = .5)






set.seed(23235)
ss <- sample(1:150, 10)
dend1 <- iris[ss, -5] %>%
  dist() %>%
  hclust("com") %>%
  as.dendrogram()
dend2 <- iris[ss, -5] %>%
  dist() %>%
  hclust("sin") %>%
  as.dendrogram()
dend12 <- dendlist(dend1, dend2)

dend12 %>% tanglegram()

dend12
dndlist
