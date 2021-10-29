library(devtools)
install_github("thackl/thacklr")
install_github("thackl/gggenomes")

fileConn<-file("results/7_ANI/2_plot/.gggenomes_done")
writeLines(c("installed"), fileConn)
close(fileConn)