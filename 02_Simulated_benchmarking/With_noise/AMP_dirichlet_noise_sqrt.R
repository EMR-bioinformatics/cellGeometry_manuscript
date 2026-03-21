####Background####

# During my AL in July 2025, Myles has made key updates to residuals and noise
# Notably, he has added graded_log_noise function and a different scaling for sqrt_noise

# need to investigate graded_log_noise and revisit sqrt noise

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/Noise_updated")

sim_sampled_dir_all <- readRDS("../simulated_dirichlet_all.rds")
load("../simulated_dirichlet.rdata")

mk <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/cellmarkers_new.rds")
mk_update <- updateMarkers(mk, nsubclass = 500, expfilter = 0.25)

####graded log noise####

grade_log_all <- list()

SD <- c(0.1, 0.25, 0.5, 0.75, 1)

for(x in SD){
  grade_log_all[[paste("SD =", x)]] <- list()
  
  for(i in paste("Rep", 1:5)){
    grade_log_all[[paste("SD =", x)]][[i]] <- graded_log_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                                         sd = x)
  }
}

#saveRDS(grade_log_all, "grade_log_noise.rds")

#####cellGeometry deconvolution#####

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

cellgeo_grade <- list()

for(i in names(grade_log_all)){
  cellgeo_grade[[i]] <- list()
  
  for(x in names(grade_log_all[[i]])){
    cellgeo_grade[[i]][[x]]<- geo_output(grade_log_all[[i]][[x]],
                                       mk_update)
  }
}

#saveRDS(cellgeo_grade, "cellgeo_grade_noise.rds")

#####DWLS#####

library(DWLS)

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Exploring_other_tools/DWLS_test/Signature.rdata")
rm(time)

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
                            colnames(propdf))] #new addition
  metric_percent <- metric_set(sim_percent, propdf)
  
  list(output = allCounts_DWLS,
       metrics_percent = metric_percent,
       time = time)
  
}

####sqrt noise####

sqrt_all <- list()

sqrt_SD <- c(10, 25, 50, 75, 100)

scale_factor <- mean(colSums(sim_sampled_dir_all$`Times 30`$`Rep 1`)) / 1e9

sqrt_SD_scale <- sqrt(sqrt_SD * scale_factor) / 2
#9.689206 15.319981 21.665724 26.534985 30.639961

for(x in sqrt_SD){
  sqrt_all[[paste("SD =", x)]] <- list()
  
  for(i in paste("Rep", 1:5)){
    sqrt_all[[paste("SD =", x)]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                                               sd = x)
  }
}

#saveRDS(sqrt_all, "sqrt_noise_correct.rds")


#####cellgeometry deconvolution#####

cellgeo_sqrt <- list()

for(i in names(sqrt_all)){
  cellgeo_sqrt[[i]] <- list()
  
  for(x in names(sqrt_all[[i]])){
    cellgeo_sqrt[[i]][[x]]<- geo_output(sqrt_all[[i]][[x]],
                                         mk_update)
  }
}

#saveRDS(cellgeo_sqrt, "cellgeo_sqrt_noise_correct.rds")

#####DWLS#####

DWLS_sqrt <- list()

DWLS_sqrt[["SD = 0"]] <- list()

for(x in names(sqrt_all)){
  DWLS_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    DWLS_sqrt[[x]][[i]] <- DWLS_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(DWLS_sqrt, "DWLS_sqrt_noise_correct.rds")

#####MuSiC#####

load("/media/gcpeac/Myles/Deconvolution/AMP_scRNAseq_data_and_annotations.RData")

celseq_counts <- as.matrix(celseq_counts)

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

library(SingleCellExperiment)
library(MuSiC)

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

music_sqrt <- list()

music_sqrt[["SD = 0"]] <- list()

for(x in names(sqrt_all)){
  music_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    music_sqrt[[x]][[i]] <- music2_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(music_sqrt, "music_sqrt_noise_correct.rds")

#####LinDeconSeq#####

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

lin_sqrt <- list()

lin_sqrt[["SD = 0"]] <- list()

for(x in names(sqrt_all)){
  lin_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    lin_sqrt[[x]][[i]] <- lin_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(lin_sqrt, "lin_sqrt_noise_correct.rds")
