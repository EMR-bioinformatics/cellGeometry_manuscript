####Background####

# Script for running deconvolution of AMP simulated data

load("/media/gcpeac/Myles/Deconvolution/AMP_scRNAseq_data_and_annotations.RData")

celseq_counts <- as.matrix(celseq_counts)

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

library(cellGeometry)
library(tictoc)
tic()
mk <- cellMarkers(celseq_counts,
                  bulkdata = tpmdata,
                  subclass = metadata$subclass,
                  cellgroup = metadata$type, dual_mean = T)
toc(log = TRUE)

mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

mk$mk_time <- mk_time
#saveRDS(mk, file = "cellmarkers_new.rds")

####generating simulated####

set.seed(3)
sim_counts <- generate_samples(mk, 25)
sim_percent <- sim_counts / rowSums(sim_counts) * 100
sim_pseudo <- simulate_bulk(mk, sim_counts)

sim_sampled_dir <- simulate_bulk(celseq_counts,
                                 sim_counts,
                                 metadata$subclass,
                                 method = "dirichlet") #times = 30 is default

# save(sim_counts, sim_percent, sim_pseudo,
#      sim_sampled, sim_sampled_dir,
#      file = "simulated_dirichlet.rdata")

sim_sampled_dir_all <- list()

sim_sampled_dir_all[["Times 30"]][["Rep 1"]] <- sim_sampled_dir

for(x in paste("Rep", 2:5)){
  sim_sampled_dir_all[["Times 30"]][[x]] <- simulate_bulk(celseq_counts,
                                         sim_counts,
                                         metadata$subclass,
                                         method = "dirichlet")
}

#saveRDS(sim_sampled_dir_all, "simulated_dirichlet_all.rds")

####cellgeometry####

cellgeo_tune_eq <- tune_deconv(mk, sim_sampled_dir,
                               sim_percent, 
                               grid = list(nsubclass = c(5, 10, 20, 25, 50, 100, 200, 500),
                                           expfilter = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3),
                                           use_filter = c(TRUE, FALSE)),
                               output = "percent",
                               metric = "RMSE",
                               arith_mean = T,
                               count_space = T, convert_bulk = F, weight_method = "equal",
                               cores = 8)

# Best tune:
# nsubclass  expfilter  use_filter  mean.Rsq
# 500        0.5        TRUE    0.9483

cellgeo_tune_eqv2 <- tune_deconv(mk, sim_sampled_dir,
                                 sim_percent, 
                                 grid = list(nsubclass = c(5, 10, 25, 50, 100, 200, 500, 800, 1000, 2000, 5000),
                                             weight_method = c("none", "equal"),
                                             expfilter = c(0.1, 0.25, 0.5, 0.75, 1)),
                                 output = "percent",
                                 metric = "RMSE",
                                 arith_mean = T,
                                 count_space = T, convert_bulk = F, 
                                 use_filter = FALSE,
                                 cores = 8)

# nsubclass  weight_method  expfilter  mean.RMSE
# 2000          equal        0.1     0.5673

mk_update <- updateMarkers(mk, nsubclass = 500, expfilter = 0.25)

geo_output <- function(times, rep){
  
  samples_sim <- sim_sampled_dir_all[[paste("Times", times)]][[paste("Rep", rep)]]
  tic()
  out <- deconvolute(mk_update, samples_sim,
                     arith_mean = T,
                     count_space = T, convert_bulk = F, 
                     use_filter = FALSE,
                     weight_method = "equal")

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
  cellgeo[["Times 30"]][[x]] <- geo_output(times = 30,
                                           rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(cellgeo, file = "cellgeo_dirichlet_output.rds")

####MuSiC1####

library(Biobase)
library(MuSiC, lib.loc = "/usr/local/lib/R/site-library")
#version 0.2.0 on workstation

sc.pheno <- data.frame(check.names = F, check.rows = F,
                       stringsAsFactors = F,
                       row.names = colnames(celseq_counts)[is.na(metadata$subclass) == FALSE],
                       SubjectName = factor(metadata$sample[is.na(metadata$subclass) == FALSE]),
                       cellType = factor(metadata$type[is.na(metadata$subclass) == FALSE],
                                         levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell')),
                       subclass = factor(metadata$subclass[is.na(metadata$subclass) == FALSE],
                                         levels = col_order))

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType",
                                         "subclass"),
                      row.names=c("SubjectName",
                                  "cellType",
                                  "subclass"))

sc.pdata <- new("AnnotatedDataFrame",
                data = sc.pheno,
                varMetadata = sc.meta)

sce <- Biobase::ExpressionSet(assayData = as.matrix(celseq_counts[ ,is.na(metadata$subclass) == FALSE]),
                              phenoData = sc.pdata)

music_out <- function(times, rep){
  samples_sim <- sim_sampled_dir_all[[paste("Times", times)]][[paste("Rep", rep)]]
  bk.pheno <- data.frame(check.names = F, check.rows = F,
                         stringsAsFactors = F,
                         row.names = colnames(samples_sim),
                         SubjectName = colnames(samples_sim))
  
  bk.meta <- data.frame(labelDescription = c("SubjectName"),
                        row.names = c("SubjectName"))
  
  bk.pdata <- new("AnnotatedDataFrame",
                  data = bk.pheno,
                  varMetadata = bk.meta)
  
  bke <- Biobase::ExpressionSet(assayData = samples_sim ,
                                phenoData = bk.pdata)
  
  tic()
  out <- music_prop(bulk.eset = bke, 
                    sc.eset = sce,
                    clusters = "subclass",
                    samples = "SubjectName")
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

music1 <- list()

for(x in paste("Rep", 1:5)){
  music1[["Times 30"]][[x]] <- music_out(times = 30,
                                         rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(music1, "Music_dirichlet_output.rds")

####MuSiC2####

detach("package:MuSiC", unload = TRUE)

#restarted R

library(SingleCellExperiment)
library(MuSiC, lib.loc =  "/home/rachel/R/x86_64-pc-linux-gnu-library/4.2")

sce <- SingleCellExperiment(list(counts = celseq_counts[ ,is.na(metadata$subclass) == FALSE]),
                            colData = metadata[is.na(metadata$subclass) == FALSE, ])

library(tictoc)
library(cellGeometry)
music2_out <- function(times, rep){
  samples_sim <- sim_sampled_dir_all[[paste("Times", times)]][[paste("Rep", rep)]]
  
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

music2 <- list()

for(x in paste("Rep", 1:5)){
  music2[["Times 30"]][[x]] <- music2_out(times = 30,
                                          rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(music2, "Music2_dirichlet_output.rds")

####DWLS####

library(DWLS)

tic()
Signature <- buildSignatureMatrixMAST(scdata = dataSC,
                                      id = labels,
                                      path = "results")
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

#save(Signature, time, file = "Signature.rdata")

DWLS_out <- function(times, rep){
  samples_sim <- sim_sampled_dir_all[[paste("Times", times)]][[paste("Rep", rep)]]
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
    if (inherits(solDWLS, 'try-error')) solDWLS <- rep(0, ncol(S))
    #try is a new addition since had problems 
    allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  }
  toc(log = TRUE) 
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  colnames(allCounts_DWLS) <- colnames(samples_sim)
  colnames(allCounts_DWLS) <- colnames(samples_sim)
  rownames(allCounts_DWLS) <- paste0("SC-", rownames(allCounts_DWLS))
  propdf <- t(allCounts_DWLS) * 100
  propdf <- propdf[ , match(colnames(sim_percent),
                            colnames(propdf))] #new addition
  metric_percent <- metric_set(sim_percent, propdf)
  
  list(output = allCounts_DWLS,
       metrics_percent = metric_percent,
       time = time)
  
}

DWLS <- list()

for(x in paste("Rep", 1:5)){
  DWLS[["Times 30"]][[x]] <- DWLS_out(times = 30,
                                      rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(DWLS, "DWLS_dirichlet_output.rds")

####LinDeconseq####

#restart R

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

library(LinDeconSeq)
library(tictoc)
tic()
markerRes <- findMarkers(dataSC,
                         phes = phes)
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
markerRes$Time <- time

#saveRDS(markerRes, "LinDeconSeq_signature.rds")

library(cellGeometry)

linDecon_out <- function(times, rep){
  samples_sim <- sim_sampled_dir_all[[paste("Times", times)]][[paste("Rep", rep)]]
  
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

Lin <- list()

for(x in paste("Rep", 1:5)){
  Lin[["Times 30"]][[x]] <- linDecon_out(times = 30,
                                         rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(Lin, "Lin_dirichlet_output.rds")









