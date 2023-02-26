library(devtools)
install_github("thackl/thacklr")
install_github("thackl/gggenomes")
install_github("eliocamp/ggnewscale@dev")

fileConn<-file("resources/CrassUS_db/.gggenomes_install_done")
writeLines(c("installed"), fileConn)
close(fileConn)
