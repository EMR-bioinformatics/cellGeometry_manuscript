load("neuron_subset.rdata") # from LinDeconSeq

library(DWLS)
library(tictoc)
tic()
Signature <- buildSignatureMatrixMAST(scdata = mat_sub,
                                      id = subcl,
                                      path = "DWLS_results_neuron")
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

#save(Signature, time, file = "DWLS_signature_neuron.rdata")

#deconvolution

load("brain_simulated_merged.rdata")

DWLS_out <- function(times, rep){
  samples_sim <- sim_sampled_merge[[paste("Times", times)]][[paste("Rep", rep)]]
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
  DWLS[["Times 1"]][[x]] <- DWLS_out(times = 1,
                                     rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(DWLS, file = "DWLS_dirichlet_output_neuron.rds")
