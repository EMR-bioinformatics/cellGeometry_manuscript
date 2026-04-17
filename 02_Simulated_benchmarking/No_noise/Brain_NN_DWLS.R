####non-neuronal####
library(zellkonverter)
library(SingleCellExperiment)

brainNN <- readH5AD("/media/lvm1/brain/99f27be8-9fac-451e-9723-9e4c7191589e.h5ad",
                    use_hdf5 = TRUE, reader = "R")

mat2 <- brainNN@assays@data$X
rownames(mat2) <- rownames(brainNN)  # need to add rownames (genes)

#dim(mat2)
# [1]  59236 888263
# 59236 * 888263 < 2^31
#FALSE

meta2 <- brainNN@colData@listData

ind_NN <- readRDS("brain_subset_index_NN.rds")

library(checkmate)
check_names(unique(as.vector(meta2$cell_type)), "strict")
# "Must have names according to R's variable naming conventions, but element 3 does not comply"

meta2$cell_type <- gsub(" ", "_", meta2$cell_type)
check_names(unique(as.vector(meta2$cell_type)), "strict")

mat_sub2 <- as.matrix(mat2[ ,ind_NN])
subcl <- meta2$cell_type[ind_NN]

colnames(mat_sub2) <- paste0("Cell", 1:length(subcl))

subcl <- factor(subcl)
donor_id <- meta2$donor_id[ind_NN]

# save(mat_sub2, donor_id, subcl,
#      file = "NN_subset.rdata")

library(DWLS)
library(tictoc)
tic()
Signature <- buildSignatureMatrixMAST(scdata = mat_sub2,
                                      id = subcl,
                                      path = "DWLS_results")
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

#save(Signature, time, file = "DWLS_signature_NN.rdata")

#deconvolution 
load("brain_simulated_merged.rdata")

library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

rownames(Signature) <- gene2symbol(rownames(Signature), ensDb_v110)

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

DWLS <- list()

for(x in paste("Rep", 1:5)){
  DWLS[["Times 1"]][[x]] <- DWLS_out(times = 1,
                                     rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(DWLS, file = "DWLS_dirichlet_output_NN.rds")