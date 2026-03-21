####Background####

# Continuation of original script 
# LinDeconSeq and MuSiC to be tested

library(cellGeometry)
sim_sampled_dir_all <- readRDS("../simulated_dirichlet_all.rds")
load("../simulated_dirichlet.rdata")

load("/media/gcpeac/Myles/Deconvolution/AMP_scRNAseq_data_and_annotations.RData")

celseq_counts <- as.matrix(celseq_counts)

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

####MuSiC####

library(SingleCellExperiment)
library(MuSiC, lib.loc =  "/home/rachel/R/x86_64-pc-linux-gnu-library/4.2")

sce <- SingleCellExperiment(list(counts = celseq_counts[ ,is.na(metadata$subclass) == FALSE]),
                            colData = metadata[is.na(metadata$subclass) == FALSE, ])

music2_out <- function(mat){
  
  samples_sim <- mat
  
  tic()
  out <- music_prop(bulk.mtx = samples_sim,
                    sc.sce = sce,
                    clusters = "subclass",
                    samples = "sample")
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  propdf <- out$Est.prop.weighted * 100
  propdf <- propdf[ , match(colnames(sim_percent), colnames(propdf))]
  
  metrics_percent <- metric_set(sim_percent, propdf)
  
  list(output = out,
       metrics_percent = metrics_percent,
       time = mk_time)
}

library(tictoc)

######shift noise#####

noise_all <- readRDS("default_noise.rds")
shift_all <- readRDS("shift_noise.rds")

music_noise <- list()

music_noise[["SD = 0"]] <- list()

for(i in names(noise_all[["shift"]])){
  music_noise[["SD = 0.5"]][[i]] <- music2_out(noise_all[["shift"]][[i]])
}

for(x in names(shift_all)){
  music_noise[[x]] <- list()
  
  for(i in names(shift_all[[x]])){
    music_noise[[x]][[i]] <- music2_out(shift_all[[x]][[i]])
  }
}

#saveRDS(music_noise, "music_shift_noise.rds")

#####add noise#####

add_all <- readRDS("add_noise.rds")

music_add <- list()

music_add[["SD = 0"]] <- list()

for(i in names(noise_all[["add"]])){
  music_add[["SD = 100"]][[i]] <- music2_out(noise_all[["add"]][[i]])
}

for(x in names(add_all)){
  music_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    music_add[[x]][[i]] <- music2_out(add_all[[x]][[i]])
  }
}

#saveRDS(music_add, "music_add_noise.rds")

music_add <- readRDS("music_add_noise.rds")

for(x in names(add_all)[5:7]){
  music_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    music_add[[x]][[i]] <- music2_out(add_all[[x]][[i]])
  }
}
#saveRDS(music_add, "music_add_noise.rds")

#####log noise#####

log_all <- readRDS("log_noise.rds")

music_log <- list()

music_log[["SD = 0"]] <- list()

for(i in names(noise_all[["log"]])){
  music_log[["SD = 0.1"]][[i]] <- music2_out(noise_all[["log"]][[i]])
}

for(x in names(log_all)){
  music_log[[x]] <- list()
  
  for(i in names(log_all[[x]])){
    music_log[[x]][[i]] <- music2_out(log_all[[x]][[i]])
  }
}

#saveRDS(music_log, "music_log_noise.rds")


#####sqrt noise#####

sqrt_all <- readRDS("sqrt_noise.rds")

music_sqrt <- list()

music_sqrt[["SD = 0"]] <- list()

for(i in names(noise_all[["sqrt"]])){
  music_sqrt[["SD = 100"]][[i]] <- music2_out(noise_all[["sqrt"]][[i]])
}

for(x in names(sqrt_all)){
  music_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    music_sqrt[[x]][[i]] <- music2_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(music_sqrt, "music_sqrt_noise.rds")

music_sqrt <- readRDS("music_sqrt_noise.rds")

for(x in names(sqrt_all)[5:8]){
  music_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    music_sqrt[[x]][[i]] <- music2_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(music_sqrt, "music_sqrt_noise.rds")

####LinDeconSeq####

#restarted R

dataSC <- celseq_counts[ , is.na(metadata$subclass) == FALSE]
meta <- metadata[is.na(metadata$subclass) == FALSE, ]

#they want a binary matrix of phenotype
library(tidyr)
phes <- as.data.frame(meta %>% pivot_wider(subclass,
                                           names_from = cell_name,
                                           values_from = cell_name))
rownames(phes) <- phes$subclass
phes <- phes[ , -1]
phes <- as.matrix(phes)
phes <- ifelse(is.na(phes), 0, 1)

markerRes <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/LinDeconSeq_signature.rds")

library(LinDeconSeq)
library(cellGeometry)
library(tictoc)

lin_out <- function(mat){
  samples_sim <- mat
  
  tic()
  fractions <- deconSeq(samples_sim , markerRes$sigMatrix$sig.mat, verbose = FALSE)
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  propdf <- fractions * 100
  colnames(propdf) <- gsub("\\.", "-", colnames(propdf))
  propdf <- propdf[ , match(colnames(sim_percent),
                            colnames(propdf))]
  
  metric_percent <- metric_set(sim_percent, propdf)
  
  list(output = fractions,
       metrics_percent = metric_percent,
       time = time)
  
}

#####shift noise####

lin_noise <- list()

lin_noise[["SD = 0"]] <- list()

for(i in names(noise_all[["shift"]])){
  lin_noise[["SD = 0.5"]][[i]] <- lin_out(noise_all[["shift"]][[i]])
}

for(x in names(shift_all)){
  lin_noise[[x]] <- list()
  
  for(i in names(shift_all[[x]])){
    lin_noise[[x]][[i]] <- lin_out(shift_all[[x]][[i]])
  }
}

#saveRDS(lin_noise, "lin_shift_noise.rds")

#####add noise#####

lin_add <- list()

lin_add[["SD = 0"]] <- list()

for(i in names(noise_all[["add"]])){
  lin_add[["SD = 100"]][[i]] <- lin_out(noise_all[["add"]][[i]])
}

for(x in names(add_all)){
  lin_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    lin_add[[x]][[i]] <- lin_out(add_all[[x]][[i]])
  }
}

#saveRDS(lin_add, "lin_add_noise.rds")

lin_add <- readRDS("lin_add_noise.rds")

for(x in names(add_all)[5:7]){
  lin_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    lin_add[[x]][[i]] <- lin_out(add_all[[x]][[i]])
  }
}

#saveRDS(lin_add, "lin_add_noise.rds")


#####log noise#####

lin_log <- list()

lin_log[["SD = 0"]] <- list()

for(i in names(noise_all[["log"]])){
  lin_log[["SD = 0.1"]][[i]] <- lin_out(noise_all[["log"]][[i]])
}

for(x in names(log_all)){
  lin_log[[x]] <- list()
  
  for(i in names(log_all[[x]])){
    lin_log[[x]][[i]] <- lin_out(log_all[[x]][[i]])
  }
}

#saveRDS(lin_log, "lin_log_noise.rds")

#####sqrt noise#####

lin_sqrt <- list()

lin_sqrt[["SD = 0"]] <- list()

for(i in names(noise_all[["sqrt"]])){
  lin_sqrt[["SD = 100"]][[i]] <- lin_out(noise_all[["sqrt"]][[i]])
}

for(x in names(sqrt_all)){
  lin_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    lin_sqrt[[x]][[i]] <- lin_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(lin_sqrt, "lin_sqrt_noise.rds")

lin_sqrt <- readRDS("lin_sqrt_noise.rds")

for(x in names(sqrt_all)[5:8]){
  lin_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    lin_sqrt[[x]][[i]] <- lin_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(lin_sqrt, "lin_sqrt_noise.rds")








