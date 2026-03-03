setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2")

load("Tabula_subset.rdata")

library(DWLS)
library(tictoc)
tic()
Signature <- buildSignatureMatrixMAST(scdata = mat_sub,
                                      id = subcl,
                                      path = "DWLS_results")
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

save(Signature, time, file = "DWLS_signature.rdata")