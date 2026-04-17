####Tabula####

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/instaprism_tabula")

#signature curation 

library(zellkonverter)
library(SingleCellExperiment)

h5 <- readH5AD("/media/lvm1/tabula/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad",
               use_hdf5 = TRUE, reader = "R")

mat <- h5@assays@data$X
rownames(mat) <- rownames(h5)

meta <- h5@colData@listData

library(InstaPrism)
library(tictoc)

tic()
ref <- refPrepare(sc_Expr = mat,
                  cell.type.labels = meta$cell_type,
                  cell.state.labels = meta$cell_type)
toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

# save(ref, mk_time,
#      file = "instaprism_tabula_ref.rdata")

# deconvolution 

library(AnnotationHub)

ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

load("../tabula_simulated_workstation2.rdata")

rownames(ref@phi.cs) <- gene2symbol(rownames(ref@phi.cs), ensDb_v110)

insta_out <- function(times, rep){
  samples_sim <- sim_sampled[[paste("Times", times)]][[paste("Rep", rep)]]
  
  tic()
  out <- InstaPrism(bulk_Expr = samples_sim, refPhi_cs = ref)
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  list(output = out,
       time = mk_time)
}

insta <- list()

for(x in paste("Rep", 1:5)){
  insta[["Times 3"]][[x]] <- insta_out(times = 3,
                                       rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(insta, file = "instaprism_tabula_output.rds")











