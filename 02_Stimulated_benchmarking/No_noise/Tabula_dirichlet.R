####Background####

# Script for running deconvolution of Tabula Sapiens simulated data 

library(zellkonverter)
library(SingleCellExperiment)

h5 <- readH5AD("/media/lvm1/tabula/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad",
               use_hdf5 = TRUE, reader = "R")

mat <- h5@assays@data$X
rownames(mat) <- rownames(h5)

meta <- h5@colData@listData

library(cellGeometry)
library(tictoc)

tic()
mk <- cellMarkers(mat, subclass = meta$cell_type, 
                    cellgroup = meta$compartment, 
                    #cores = 8,
                    dual_mean = TRUE)

toc(log = TRUE)
#598.683 sec for 8 cores

mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

mk$mk_time <- mk_time
mk$mk_time_8core <- 598.683

library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

mk <- gene2symbol(mk, ensDb_v110)

# saveRDS(mk, file = "tabula_markers_dualmeans_rerun.rds")

####Generate simulated data####

library(DelayedArray)
setAutoBlockSize(2e9)

set.seed(3)
sim_counts <- generate_samples(mk, 25)
sim_percent <- sim_counts / rowSums(sim_counts) * 100
sim_pseudo <- simulate_bulk(mk, sim_counts)

sim_sampled <- list()

sim_sampled[["Times 3"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampled[["Times 3"]][[x]] <- simulate_bulk(mat, 
                                                 sim_counts,
                                                 meta$cell_type,
                                                 times = 3,
                                                 method = "dirichlet")
}

for(x in paste("Rep", 1:5)){
  rownames(sim_sampled[["Times 3"]][[x]]) <- gene2symbol(rownames(sim_sampled[["Times 3"]][[x]]), ensDb_v110)
}

# save(sim_counts, sim_percent, sim_pseudo, sim_sampled,
#      file = "tabula_simulated_workstation2.rdata")

####CellGeometry####

cellgeo_tune_eq <- tune_deconv(mkv2,  sim_sampled$`Times 3`$`Rep 1`,
                               sim_percent, 
                               grid = list(nsubclass = c(5, 10, 20, 25, 50, 100, 200, 500),
                                           expfilter = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2),
                                           use_filter = c(TRUE, FALSE)),
                               output = "percent",
                               metric = "RMSE",
                               arith_mean = T,
                               count_space = T, convert_bulk = F, weight_method = "equal",
                               cores = 8)
# Best tune:
#   nsubclass  expfilter  use_filter  mean.RMSE
# 500        0.5       FALSE    0.04986

mk_update <- updateMarkers(mk, nsubclass = 500, expfilter = 0.5)

geo_output <- function(times, rep){
  
  samples_sim <- sim_sampled[[paste("Times", times)]][[paste("Rep", rep)]]
  tic()
  out <- deconvolute(mk_update, samples_sim,
                     arith_mean = T,
                     count_space = T, convert_bulk = F, 
                     use_filter = FALSE, 
                     #cores = 1,
                     cores = 8,
                     weight_method = "equal")
  #original use_filter = T
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  metrics_counts <- metric_set(sim_counts, out$subclass$output/times)
  metrics_percent <- metric_set(sim_percent, out$subclass$percent)
  
  list(output = out,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = mk_time)
  
}

cellgeo <- list()

for(x in paste("Rep", 1:5)){
  cellgeo[["Times 3"]][[x]] <- geo_output(times = 3,
                                          rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(cellgeo, file = "cellgeo_dirichlet_output_500.rds") # 1 core
# saveRDS(cellgeo, file = "cellgeo_dirichlet_output_500_8cores.rds")

tic()
mk_update <- updateMarkers(mk, nsubclass = 200, expfilter = 0.5)
toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
mk_update$mk_time <- mk_time

# saveRDS(mk_update, "tabula_markers_dualmeans_rerun_update.rds")

cellgeo <- list()

for(x in paste("Rep", 1:5)){
  cellgeo[["Times 3"]][[x]] <- geo_output(times = 3,
                                          rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(cellgeo, file = "cellgeo_dirichlet_output_200_8cores_rerun.rds")

####MuSiC####

library(MuSiC)
library(SingleCellExperiment)

library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

rownames(mat) <- gene2symbol(rownames(mat), ensDb_v110)

sce <- SingleCellExperiment(list(counts = mat[ , is.na(meta$cell_type) == FALSE]),
                            colData = data.frame("subclass" = meta$cell_type[is.na(meta$cell_type) == FALSE],
                                                 "donor_id" = meta$donor_id[is.na(meta$cell_type) == FALSE]))

music2_out <- function(times, rep){
  samples_sim <- sim_sampled[[paste("Times", times)]][[paste("Rep", rep)]]
  
  tic()
  out <- music_prop(bulk.mtx = samples_sim,
                    sc.sce = sce,
                    clusters = "subclass",
                    samples = "donor_id")
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  propdf <- out$Est.prop.weighted * 100
  
  common <- intersect(colnames(sim_percent), colnames(propdf))
  
  propdf <- propdf[ , colnames(propdf) %in% common]
  
  missing <- setdiff(colnames(sim_percent), colnames(propdf))
  
  propdf <- as.data.frame(propdf)
  
  if(length(missing) > 0){
    for(i in missing){
      propdf[ , i] <- 0
    }
  }
  
  propdf <- propdf[ , match(colnames(sim_percent), colnames(propdf))]
  
  metrics_percent <- metric_set(sim_percent, propdf)
  
  list(output = out,
       metrics_percent = metrics_percent,
       time = mk_time)
}

music <- list()

for(x in paste("Rep", 1:5)){
  music[["Times 3"]][[x]] <- music2_out(times = 3,
                                        rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(music, file = "music_dirichlet_output.rds")






