####non-neuronal####
library(zellkonverter)
library(SingleCellExperiment)

brainNN <- readH5AD("/media/lvm1/brain/99f27be8-9fac-451e-9723-9e4c7191589e.h5ad",
                    use_hdf5 = TRUE, reader = "R")

mat2 <- brainNN@assays@data$X
rownames(mat2) <- rownames(brainNN)  

#dim(mat2)
# [1]  59236 888263
# 59236 * 888263 < 2^31
#FALSE

meta2 <- brainNN@colData@listData

slim <- function(subclass, maxlen = 200) {
  tab <- table(subclass)
  ok <- tab >= 10 & tab <= maxlen
  sam <- which(tab > maxlen)
  ind <- which(subclass %in% levels(subclass)[ok])
  if (length(sam) > 0) {
    ind2 <- unlist(lapply(sam, function(i) {
      subcl <- which(subclass == levels(subclass)[i])
      sample(subcl, maxlen)
    }))
    ind <- sort(c(ind, ind2))
  }
  ind
}

ind_NN <- slim(meta2$cell_type, 100) 

#saveRDS(ind_NN, "brain_subset_index_NN.rds")

library(checkmate)
check_names(unique(as.vector(meta2$cell_type)), "strict")
# "Must have names according to R's variable naming conventions, but element 3 does not comply"

meta2$cell_type <- gsub(" ", "_", meta2$cell_type)
check_names(unique(as.vector(meta2$cell_type)), "strict")
#TRUE

meta_sub2 <- data.frame("Cell_type" = meta2$cell_type [ind_NN],
                       "Cell" = meta2$observation_joinid[ind_NN])

meta_sub2$Cell <- paste0("Cell", 1:nrow(meta_sub2))

library(tidyverse)
phes2 <- as.data.frame(meta_sub2 %>% pivot_wider(Cell_type,
                                               names_from = Cell,
                                               values_from = Cell))
rownames(phes2) <- phes2$Cell_type
phes2 <- phes2[ , -1]
phes2 <- as.matrix(phes2)
phes2 <- ifelse(is.na(phes2), 0, 1)

#table(colSums(phes2))
#all 1s so fine

mat_sub2 <- as.matrix(mat2[ ,ind_NN])

colnames(mat_sub2) <- paste0("Cell", 1:nrow(meta_sub2))

library(LinDeconSeq)
library(tictoc)
tic()
markerRes <- findMarkers(mat_sub2,
                         phes = phes2)
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
markerRes$Time <- time

# [INFO] 1100 samples and 59236 genes in the reference profile
# [INFO] 3024 candidate cell type-specific genes detected with q.cut 0.01
# [INFO] Select seed genes and random permutating...
# [INFO] Assigning cell type-specific genes to each cell subset
# [INFO] 60 low confidence marker genes have been filtered out.                                                                                                            
# [INFO] Optimizing cell type-specific genes to derive signature matrix...
# [INFO] Group size -> 88, condition number -> 3.875581

#saveRDS(markerRes, file = "LinDeconSeq_sig_NN.rds")

# deconvolution 

sig <- sig_obj$sigMatrix$sig.mat
library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]
rownames(sig) <- gene2symbol(rownames(sig), ensDb_v110)

linDecon_out <- function(times, rep){
  samples_sim <- sim_sampled_merge[[paste("Times", times)]][[paste("Rep", rep)]]
  
  tic()
  fractions <- deconSeq(samples_sim , sig, verbose = FALSE)
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  list(output = fractions,
       time = time)
  
}

lin <- list()

for(x in paste("Rep", 1:5)){
  lin[["Times 1"]][[x]] <- linDecon_out(times = 1,
                                        rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(lin, file = "Lin_dirichlet_output_NN.rds")


