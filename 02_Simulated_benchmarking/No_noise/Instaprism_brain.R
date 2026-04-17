####neuron####

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/Instaprism_brain")

#signature curation 

library(zellkonverter)
library(SingleCellExperiment)

brain <- readH5AD("/media/lvm1/brain/c2f66cd5-4ff4-4578-876c-55783a57cf8f.h5ad",
                  use_hdf5 = TRUE, reader = "R")

mat <- brain@assays@data$X
rownames(mat) <- rownames(brain)  
meta <- brain@colData@listData

library(InstaPrism)
library(tictoc)

tic()
ref <- refPrepare(sc_Expr = mat,
                  cell.type.labels = meta$roi,
                  cell.state.labels = meta$roi)
toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

# save(ref, mk_time,
#      file = "instaprism_neuron_ref.rdata")

library(AnnotationHub)
library(cellGeometry)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

load("../brain_simulated_merged.rdata")

rownames(ref@phi.cs) <- gene2symbol(rownames(ref@phi.cs), ensDb_v110)

insta_out <- function(times, rep){
  samples_sim <- sim_sampled_merge[[paste("Times", times)]][[paste("Rep", rep)]]
  
  tic()
  out <- InstaPrism(bulk_Expr = samples_sim, refPhi_cs = ref)
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  list(output = out,
       time = mk_time)
}

insta <- list()

for(x in paste("Rep", 1:3)){
  insta[["Times 1"]][[x]] <- insta_out(times = 1,
                                       rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(insta, file = "instaprism_neuron_output.rds")

####non-neuron####

brainNN <- readH5AD("/media/lvm1/brain/99f27be8-9fac-451e-9723-9e4c7191589e.h5ad",
                    use_hdf5 = TRUE, reader = "R")

mat2 <- brainNN@assays@data$X
rownames(mat2) <- rownames(brainNN)

meta2 <- brainNN@colData@listData

tic()
ref_NN <- refPrepare(sc_Expr = mat2,
                     cell.type.labels = meta2$cell_type,
                     cell.state.labels = meta2$cell_type)
toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

# save(ref_NN, mk_time,
#      file = "instaprism_NN_ref.rdata")

rownames(ref_NN@phi.cs) <- gene2symbol(rownames(ref_NN@phi.cs), ensDb_v110)

insta_out <- function(times, rep){
  samples_sim <- sim_sampled_merge[[paste("Times", times)]][[paste("Rep", rep)]]
  
  tic()
  out <- InstaPrism(bulk_Expr = samples_sim, refPhi_cs = ref_NN)
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  list(output = out,
       time = mk_time)
}

insta <- list()

for(x in paste("Rep", 1:5)){
  insta[["Times 1"]][[x]] <- insta_out(times = 1,
                                       rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(insta, file = "instaprism_NN_output.rds")

