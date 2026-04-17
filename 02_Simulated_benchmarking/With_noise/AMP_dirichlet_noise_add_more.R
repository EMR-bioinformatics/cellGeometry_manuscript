music_add <- readRDS("music_add_noise.rds")
cellgeo_add <- readRDS("cellgeo_add_noise.rds")
DWLS_add <- readRDS("DWLS_add_noise.rds")
Lin_add <- readRDS("lin_add_noise.rds")
add_all <- readRDS("add_noise.rds")
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

####add SD more####

add_all[["SD = 2000"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 2000"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 2000)
}

add_all[["SD = 3000"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 3000"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 3000)
}

add_all[["SD = 4000"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 4000"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 4000)
}

add_all[["SD = 5000"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 5000"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                            sd = 5000)
}

#saveRDS(add_all, file = "add_noise.rds")

####cellGeometry####

geo_output <- function(mat, mk){
  
  samples_sim <- mat
  tic()
  out <- deconvolute(mk, samples_sim,
                     arith_mean = T,
                     count_space = T, convert_bulk = F, 
                     use_filter = FALSE,
                     weight_method = "equal")
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  metrics_counts <- metric_set(sim_counts, out$subclass$output/30)
  metrics_percent <- metric_set(sim_percent, out$subclass$percent)
  
  list(output = out,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = mk_time)
  
}

mk <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/cellmarkers_new.rds")
mk_update <- updateMarkers(mk, nsubclass = 1000, expfilter = 2.5)

for(x in paste0("SD = ", seq(1000, 5000, 1000))){
  cellgeo_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    cellgeo_add[[x]][[i]] <- geo_output(add_all[[x]][[i]], mk_update)
  }
}

#saveRDS(cellgeo_add, file = "cellgeo_add_noisev2.rds")

#####music####

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

library(MuSiC)
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts = celseq_counts[ ,is.na(metadata$subclass) == FALSE]),
                                 colData = metadata[is.na(metadata$subclass) == FALSE, ])

for(x in paste0("SD = ", seq(1000, 5000, 1000))){
  music_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    music_add[[x]][[i]] <- music2_out(add_all[[x]][[i]])
  }
}

#saveRDS(music_add, file = "music_add_noisev2.rds")

####DWLS####

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Exploring_other_tools/DWLS_test/Signature.rdata")
rm(time)
library(DWLS)

DWLS_out <- function(mat){
  samples_sim <- mat
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
  colnames(allCounts_DWLS) <- colnames(samples_sim)
  rownames(allCounts_DWLS) <- paste0("SC-", rownames(allCounts_DWLS))
  propdf <- t(allCounts_DWLS) * 100
  propdf[is.na(propdf) == TRUE] <- 0
  propdf <- propdf[ , match(colnames(sim_percent),
                            colnames(propdf))] 
  metric_percent <- metric_set(sim_percent, propdf)
  
  list(output = allCounts_DWLS,
       metrics_percent = metric_percent,
       time = time)
  
}

for(x in paste0("SD = ", seq(1000, 5000, 1000))){
  DWLS_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    DWLS_add[[x]][[i]] <- DWLS_out(add_all[[x]][[i]])
  }
}

#saveRDS(DWLS_add, file = "DWLS_add_noisev2.rds")

# LinDeconSeq

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
  
  metric_percent <- metric_set(sim_percent, as.matrix(propdf))
  
  list(output = fractions,
       metrics_percent = metric_percent,
       time = time)
  
}

for(x in paste0("SD = ", seq(1000, 5000, 1000))){
  Lin_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    Lin_add[[x]][[i]] <- lin_out(add_all[[x]][[i]])
  }
}

#saveRDS(Lin_add, file = "Lin_add_noisev2.rds")

#combine

cellgeo_add <- readRDS("cellgeo_add_noisev2.rds")
DWLS_add <- readRDS("DWLS_add_noisev2.rds")
music_add <- readRDS("music_add_noisev2.rds")
Lin_add <- readRDS("Lin_add_noisev2.rds")

cellgeo_baseline <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/cellgeo_dirichlet_output500.rds")
DWLS_baseline <- readRDS("../DWLS_dirichlet_output.rds")
music_baseline <- readRDS("../music_dirichlet_output.rds")
lin_baseline <- readRDS("../Lin_dirichlet_output.rds")

percent_retrieve <- function(list, method){
  df <- data.frame()
  for(i in names(list)){
    for(x in names(list[[i]])){
      temp <- as.data.frame(list[[i]][[x]]$metrics_percent)
      temp <- temp[ , c("Rsq", "RMSE")]
      temp$Rep <- as.numeric(gsub("Rep ", "", x))
      temp$Noise <- as.numeric(gsub("SD = ", "", i))
      temp <- add_column(temp,
                         "Cluster" = rownames(temp),
                         .before = 1)
      df <- rbind(df,
                  temp)
    }
  }
  df$Method <- method
  
  df
}

add_noise <- rbind(percent_retrieve(cellgeo_add, "CellGeometry"),
                   percent_retrieve(music_add, "MuSiC"),
                   percent_retrieve(DWLS_add, "DWLS"),
                   percent_retrieve(Lin_add, "LinDeconSeq"))

add_noise <- add_noise[add_noise$Noise != 0, ]

add_noise <- rbind(add_noise,
                   percent_retrieve(cellgeo_baseline, "CellGeometry"),
                   percent_retrieve(music_baseline, "MuSiC"),
                   percent_retrieve(DWLS_baseline, "DWLS"),
                   percent_retrieve(lin_baseline, "LinDeconSeq"))

add_noise$Noise[is.na(add_noise$Noise) == TRUE] <- 0

add_noise$Type <- "Gaussian noise"

add_noise$Method <- factor(add_noise$Method,
                           levels = c("CellGeometry",
                                      "CellGeometry 2-npass",
                                      "MuSiC",
                                      "DWLS",
                                      "LinDeconSeq"))

#saveRDS(add_noise , file = "add_noise_results.rds")

