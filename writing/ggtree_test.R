library(ggtree)
library(ape)
library(dplyr)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(phytools)
library(ggpubr)
library(cowplot)
library(ggstance)

setwd("~/encode_home/projects/crAssUS/")

# read tree
t <- read.tree("writing/repr_portal.nwk")
t$tip.label <- sapply(strsplit(t$tip.label,"|", fixed=T), `[`, 1)

# read annot
annot_file <- read.csv('writing/repr_portal_annot.txt', sep = '\t', header = T)


# group families to color them later
dat4 <- annot_file %>% select(c("genome", "family"))
dat4 <- aggregate(.~family, dat4, FUN=paste, collapse=",")
clades <- lapply(dat4$genome, function(x){unlist(strsplit(x,split=","))})
names(clades) <- dat4$family
t1 <- groupOTU(t, clades, 'Family')



p <- ggtree(tr=t1, aes(color=Family), layout = "fan", size=0.1) + 
  scale_color_manual(values = c ("Jelitoviridae" = "#2ca02c",
                                 "Crevaviridae" = "#9467bd",
                                 "Intestiviridae" = "#ff7f0e",
                                 "Steigviridae" = "#1f77b4",
                                 "Suoliviridae" = "#d62728", 
                                 "Tinaiviridae" = "#8c564b")) #+
  # geom_tiplab(aes(subset=(grepl('Howler',label,fixed=TRUE)==TRUE)),
  #             size=0.8)


p1 <- p +
  geom_fruit(
    data=annot_file,
    geom=geom_star,
    mapping=aes(y=genome, fill=sample_type, starshape=shape),
    size=0.3,
    colour = "white",
    starstroke=0,
    #pwidth=0.1,
    inherit.aes = FALSE,
    # grid.params=list(
    #   linetype=0,
    #   size=0.2,
    # )
    
  ) +
  scale_fill_manual(
    name="Source",
    values=c("ancient"="magenta", "modern"="orange", "iceman"="lightseagreen","NH_primates"="saddlebrown"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=4,
                       override.aes=list(starshape=c("ancient"=15,
                                                     "modern"=15,
                                                     "iceman"=15,
                                                     "NH_primates"=15),
                                         size=2)
    ),
    na.translate=FALSE,
  ) +
  scale_starshape_manual(
    values=c("circle"=15),
    guide="none"
  ) +
  new_scale_fill()

p2 <- p1 +
  geom_fruit(
    data=annot_file,
    geom=geom_tile,
    mapping=aes(y=genome, fill=coding),
    width=0.1,
    color="white",
    pwidth=0.1,
    offset=0.05
  ) +
  scale_fill_manual(
    name="Coding",
    values=c("grey80", "orchid4", "chartreuse4"),
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=3
    )
  ) +
  new_scale_fill()

portal <- p2 +
  geom_fruit(data=annot_file, geom=geom_bar,
             mapping=aes(y=genome, x=compl, fill=final_family),
             pwidth=0.2, 
             orientation="y", 
             stat="identity",
             grid.params=list(
               linetype=2,
               size=0.2,
             )
) + scale_fill_manual(values = c ("Jelitoviridae" = "#2ca02c",
                                   "Crevaviridae" = "#9467bd",
                                   "Intestiviridae" = "#ff7f0e",
                                   "Steigviridae" = "#1f77b4",
                                   "Suoliviridae" = "#d62728", 
                                   "Tinaiviridae" = "#8c564b"),
                      guide="none")

terl <- terl +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0.02, "cm"),
        plot.margin=margin(0, 0, 0, 0, "mm"))

# mcp
# ggsave("writing/ggtree/repr_mcp", device="tiff", dpi=300)


prow <- plot_grid(
  terl + theme(legend.position="none"),
  mcp + theme(legend.position="none"),
  portal + theme(legend.position="none"),
  NULL,
  align = 'vh',
  labels = c("TerL", "MCP", "Portal", ""),
  hjust = -1,
  nrow = 1,
  rel_widths = c(1, 1, 1, .3)
)

# https://wilkelab.org/cowplot/articles/shared_legends.html
# https://stackoverflow.com/questions/41570188/plot-legend-in-an-empty-panel-with-cowplot-ggplot2
legend <- get_legend(
                    terl + theme(legend.box.margin = margin(0, 0, 0, 0))

                    )
# legend <- get_legend(
#                     terl + theme(legend.box.margin = margin(0, 0, 0, 0),
#                                  legend.position = "bottom")
#                     )

pdf("writing/ggtree/Fig_2.pdf", family = "ArialMT", width = 10, height = 5, pointsize = 18)
prow + draw_grob(legend, 3/3.3, 0, .1/3.3, 1)  # legend, 3/3.3, 0, .3/3.3, .3
#plot_grid(prow, legend, ncol=1)
invisible(dev.off())

dev.off()
#ggsave("writing/ggtree/repr_all", device="svg", dpi=400, width = 2200, height = 800, units = "px")
# 

###############
##### GPD #####
###############

# read tree
t <- read.tree("writing/repr_gpd_TerL.nwk")
t$tip.label <- sapply(strsplit(t$tip.label,"|", fixed=T), `[`, 1)

# read annot
annot_file <- read.csv('writing/repr_gpd_TerL_annot.txt', sep = '\t', header = T)


# group families to color them later
dat4 <- annot_file %>% select(c("genome", "family"))
dat4 <- aggregate(.~family, dat4, FUN=paste, collapse=",")
clades <- lapply(dat4$genome, function(x){unlist(strsplit(x,split=","))})
names(clades) <- dat4$family
t1 <- groupOTU(t, clades, 'Family')


# read counts data
counts <- read.csv('writing/gpd_source_counts_repr.txt', sep = '\t', header = T)
sum_value <- counts %>%
  group_by(genome) %>%
  summarize(total = sum(count))
sum_value[sum_value == 0] <- NA


#read tippoints file with the reference crassphages
tippoints <- read.csv('writing/gpd_tippoints_reference.txt', sep = '\t', header = T)


p <- ggtree(tr=t1, aes(color=Family), size=0.1, show.legend=F)  + 
  scale_color_manual(values = c ("Jelitoviridae" = "#2ca02c",
                                 "Crevaviridae" = "#9467bd",
                                 "Intestiviridae" = "#ff7f0e",
                                 "Steigviridae" = "#1f77b4",
                                 "Suoliviridae" = "#d62728", 
                                 "Tinaiviridae" = "#8c564b")) +
  new_scale_color()

p1 <- p %<+% tippoints + geom_tippoint(aes(color=family, shape=family), size=0.1, show.legend = F) +
  scale_color_manual(values = c ("Jelitoviridae" = "#2ca02c",
                                 "Crevaviridae" = "#9467bd",
                                 "Intestiviridae" = "#ff7f0e",
                                 "Steigviridae" = "#1f77b4",
                                 "Suoliviridae" = "#d62728", 
                                 "Tinaiviridae" = "#8c564b"),
                      na.translate=FALSE, 
                      ) + scale_shape_manual(values = c ("Jelitoviridae" = 20,
                                                         "Crevaviridae" = 20,
                                                         "Intestiviridae" = 20,
                                                         "Steigviridae" = 20,
                                                         "Suoliviridae" = 20,
                                                         "Tinaiviridae" = 20))

p3 <- facet_plot(p1, panel = 'CrAssUS and GPD agreement', data = counts, 
                 geom = geom_barh, 
                 mapping = aes(x = count, fill = as.factor(source)), 
                 stat='identity',
                 show.legend=F) 
p4 <- facet_plot(p3, panel = 'CrAssUS and GPD agreement', data = sum_value, 
                 geom = geom_text,
                 mapping=aes(x=total+20, label=round(total)),size=1.2, color="gray30",
                 show.legend=F) 


pdf("writing/ggtree/Fig_4.pdf", family = "ArialMT", width = 15 / 2.54, height = 28/ 2.54, pointsize = 4)
#plot_grid(prow, legend, rel_widths = c(1, 1, 1, .2))
#prow + draw_grob(legend, 3.1/3.3, .3/1, .3/3.3, .3, scale = 0.5)  # legend, 3/3.3, 0, .3/3.3, .3
p4
invisible(dev.off())
