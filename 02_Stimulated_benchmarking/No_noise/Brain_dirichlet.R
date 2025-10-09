####Background####

# Script for running deconvolution of brain simulated data

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain")

library(zellkonverter)
library(SingleCellExperiment)
library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

brain <- readH5AD("/media/lvm1/brain/c2f66cd5-4ff4-4578-876c-55783a57cf8f.h5ad",
                  use_hdf5 = TRUE, reader = "R")

mat <- brain@assays@data$X
rownames(mat) <- rownames(brain)  # need to add rownames (genes)
meta <- brain@colData@listData

library(cellGeometry)
library(tictoc)
tic()
mk <- cellMarkers(mat, subclass = meta$roi,
                  cellgroup = meta$supercluster_term, 
                  dual_mean = TRUE, cores = 8)
toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

mk$mk_time <- mk_time

mk <- gene2symbol(mk, ensDb_v110)

#saveRDS(mk, "brain_cellmarkers.rds")

# non-neuronal cells
brainNN <- readH5AD("/media/lvm1/brain/99f27be8-9fac-451e-9723-9e4c7191589e.h5ad",
                    use_hdf5 = TRUE, reader = "R")

mat2 <- brainNN@assays@data$X
rownames(mat2) <- rownames(brainNN)  # need to add rownames (genes)
meta2 <- brainNN@colData@listData

tic()
mkNN <- cellMarkers(mat2, subclass = meta2$cell_type, #note roi some overlap with neuron but this is non-neuron cell types
                    cellgroup = meta2$supercluster_term,
                    dual_mean = TRUE, cores = 8)
toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

mkNN$mk_time <- mk_time
mkNN <- gene2symbol(mkNN, ensDb_v110)
#saveRDS(mkNN, "brain_cellmarkers_NN.rds")

#generate simulated
set.seed(3)
sim_counts <- generate_samples(mk, 30)
sim_percent <- sim_counts / rowSums(sim_counts) * 100

sim_sampled <- list()

sim_sampled[["Times 1"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampled[["Times 1"]][[x]] <- simulate_bulk(mat, 
                                                 sim_counts,
                                                 meta$roi,
                                                 times = 1,
                                                 method = "dirichlet")
}

for(x in paste("Rep", 1:5)){
  rownames(sim_sampled[["Times 1"]][[x]]) <- gene2symbol(rownames(sim_sampled[["Times 1"]][[x]]), ensDb_v110)
}

# save(sim_counts, sim_percent, sim_sampled,
#      file = "brain_simulated.rdata")

#do the same with non-neuronal

set.seed(3)
sim_countsNN <- generate_samples(mkNN, 30)
sim_percentNN <- sim_countsNN / rowSums(sim_countsNN) * 100

sim_sampledNN <- list()

sim_sampledNN[["Times 1"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampledNN[["Times 1"]][[x]] <- simulate_bulk(mat2, 
                                                   sim_countsNN,
                                                   meta2$cell_type,
                                                   times = 1,
                                                   method = "dirichlet")
}

for(x in paste("Rep", 1:5)){
  rownames(sim_sampledNN[["Times 1"]][[x]]) <- gene2symbol(rownames(sim_sampledNN[["Times 1"]][[x]]), ensDb_v110)
}

# save(sim_countsNN, sim_percentNN, sim_sampledNN,
#      file = "brain_simulated_NN.rdata")

# merge

tic()
mkm <- mergeMarkers(mk, mkNN, transform = "none")
toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
mkm$mk_time <- mk_time

#saveRDS(mkm, "brain_cellmarkers_merged.rds")

identical(rownames(sim_sampled$`Times 1`$`Rep 1`), rownames(sim_sampledNN$`Times 1`$`Rep 1`))
#TRUE

sim_counts_merge <- cbind(sim_counts, sim_countsNN)
sim_percent_merge <- sim_counts_merge / rowSums(sim_counts_merge) * 100

sim_sampled_merge <- list()

sim_sampled_merge[["Times 1"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampled_merge[["Times 1"]][[x]] <- sim_sampled[["Times 1"]][[x]] + sim_sampledNN[["Times 1"]][[x]]
}

# save(sim_counts_merge, sim_percent_merge, sim_sampled_merge,
#      file = "brain_simulated_merged.rdata")

####CellGeometry####

res <- tune_deconv(mkm, sim_sampled_merge$`Times 1`$`Rep 1`, sim_counts_merge,
                   metric = "RMSE",
                   grid = list(nsubclass = c(10, 20, 50, 100, 200, 500),
                               expfilter = c(0, 0.05, 0.1, 0.2, 0.5, 1),
                               weight_method = c("none", "equal")),
                   arith_mean = T, count_space = T, convert_bulk = F,
                   use_filter = F, cores = 8)

# Best tune:
# nsubclass  expfilter  weight_method  mean.RMSE
# 500        0.2          equal      918.8

plot_tune(res, xvar = "nsubclass", metric = "RMSE")

mk_update <- updateMarkers(mkm, nsubclass = 500, expfilter = 0.2) 

geo_output <- function(times, rep){
  
  samples_sim <- sim_sampled_merge[[paste("Times", times)]][[paste("Rep", rep)]]
  tic()
  out <- deconvolute(mk_update, samples_sim,
                     arith_mean = T,
                     count_space = T, convert_bulk = F, 
                     use_filter = FALSE,
                     weight_method = "equal",
                     cores = 8)
  #original use_filter = T
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  metrics_counts <- metric_set(sim_counts_merge, out$subclass$output)
  metrics_percent <- metric_set(sim_percent_merge, out$subclass$percent)
  
  list(output = out,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = mk_time)
  
}

cellgeo <- list()

for(x in paste("Rep", 1:5)){
  cellgeo[["Times 1"]][[x]] <- geo_output(times = 1,
                                          rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(cellgeo, file = "cellgeo_dirichlet_output.rds")

tic()
mk_updatev2 <- updateMarkers(mkm, nsubclass = 200, expfilter = 0.2) 
toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
mk_updatev2$mk_time_update <- mk_time

#saveRDS(mk_updatev2, "brain_cellmarkers_merged_update.rds")

cellgeo <- list()

for(x in paste("Rep", 1:5)){
  cellgeo[["Times 1"]][[x]] <- geo_output(times = 1,
                                          rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(cellgeo, file = "cellgeo_dirichlet_output_rerun.rds")

####MuSiC####

library(MuSiC)

#non-neuronal first
sce <- SingleCellExperiment(list(counts = mat2),
                            colData = data.frame("subclass" = meta2$cell_type,
                                                 "donor_id" = meta2$donor_id))

music2_out <- function(times, rep){
  samples_sim <- sim_sampled_merge[[paste("Times", times)]][[paste("Rep", rep)]]
  
  tic()
  out <- music_prop(bulk.mtx = samples_sim,
                    sc.sce = sce,
                    clusters = "subclass",
                    samples = "donor_id")
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  # propdf <- out$Est.prop.weighted * 100
  # 
  # common <- intersect(colnames(sim_percentNN), colnames(propdf))
  # 
  # propdf <- propdf[ , colnames(propdf) %in% common]
  # 
  # missing <- setdiff(colnames(sim_percentNN), colnames(propdf))
  # 
  # propdf <- as.data.frame(propdf)
  # 
  # if(length(missing) > 0){
  #   for(i in missing){
  #     propdf[ , i] <- 0
  #   }
  # }
  # 
  # propdf <- propdf[ , match(colnames(sim_percentNN), colnames(propdf))]
  # 
  # metrics_percent <- metric_set(sim_percentNN, propdf)
  
  list(output = out,
       #metrics_percent = metrics_percent, #to avoid any errors - to fill later
       time = mk_time)
}

musicNN <- list()


musicNN[["Times 1"]][["Rep 1"]] <- music2_out(times = 1,
                                              rep = 1)

#saveRDS(musicNN, file = "music_dirichlet_output_NN.rds")

for(x in paste("Rep", 2:5)){
  musicNN[["Times 1"]][[x]] <- music2_out(times = 1,
                                          rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(musicNN, file = "music_dirichlet_output_NN.rds")


#when running neuronal - also make sure to run individually

sce_neuron <- SingleCellExperiment(list(counts = mat),
                                   colData = data.frame("subclass" = meta$roi,
                                                        "donor_id" = meta$donor_id))

music2_out_neuron <- function(times, rep){
  samples_sim <- sim_sampled_merge[[paste("Times", times)]][[paste("Rep", rep)]]
  
  tic()
  out <- music_prop(bulk.mtx = samples_sim,
                    sc.sce = sce_neuron,
                    clusters = "subclass",
                    samples = "donor_id")
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  list(output = out,
       time = mk_time)
}


music_neuron <- list()

music_neuron[["Times 1"]][["Rep 1"]] <- music2_out_neuron(times = 1,
                                                          rep = 1)
#saveRDS(music_neuron, file = "music_dirichlet_output.rds")

music_neuron[["Times 1"]][["Rep 2"]] <- music2_out_neuron(times = 1,
                                                          rep = 2)
#saveRDS(music_neuron, file = "music_dirichlet_output.rds")

music_neuron[["Times 1"]][["Rep 3"]] <- music2_out_neuron(times = 1,
                                                          rep = 3)
#saveRDS(music_neuron, file = "music_dirichlet_output.rds")


