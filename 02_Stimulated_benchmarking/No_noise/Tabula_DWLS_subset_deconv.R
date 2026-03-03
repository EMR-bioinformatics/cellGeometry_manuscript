setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2")

load("tabula_simulated_workstation2.rdata")

#load("DWLS_signature.rdata")

# library(AnnotationHub)
# ah <- AnnotationHub()
# ensDb_v110 <- ah[["AH113665"]]
# library(cellGeometry)
# rownames(Signature) <- gene2symbol(rownames(Signature), ensDb_v110)

#save(Signature, time, file = "DWLS_signature_converted.rdata")

load("DWLS_signature_converted.rdata")

DWLS_out <- function(times, rep){
  samples_sim <- sim_sampled[[paste("Times", times)]][[paste("Rep", rep)]]
  allCounts_DWLS<-NULL
  
  tic() 
  for(j in 1:(dim(samples_sim)[2])){
    S <- Signature
    Bulk<-samples_sim[, j]
    names(Bulk)<- rownames(samples_sim)
    Genes <- intersect(rownames(S), names(Bulk))
    B <- Bulk[Genes]
    S <- S[Genes,]
    solDWLS<-try(solveDampenedWLS(S,B), silent = TRUE)
    if (inherits(solDWLS, 'try-error')) solDWLS <- rep(NA, ncol(S))
    #try is a new addition since had problems with times = 1, rep 4
    allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  }
  toc(log = TRUE) 
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  colnames(allCounts_DWLS) <- colnames(samples_sim)

  list(output = allCounts_DWLS,
       time = time)
  
}

library(DWLS)
library(tictoc)

DWLS <- list()

for(x in paste("Rep", 1:5)){
  DWLS[["Times 3"]][[x]] <- DWLS_out(times = 3,
                                     rep = as.numeric(gsub("Rep ", "", x)))
}

saveRDS(DWLS, file = "DWLS_dirichlet_output.rds")




