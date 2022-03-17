library(devtools)
install_github("thackl/thacklr")
install_github("thackl/gggenomes")

fileConn<-file("resources/crassus_dependencies/.gggenomes_install_done")
writeLines(c("installed"), fileConn)
close(fileConn)
