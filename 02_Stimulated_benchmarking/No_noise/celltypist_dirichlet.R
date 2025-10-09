####Background####

# Script for running deconvolution of cell typist blood simulated data 

typist <- readRDS("/media/lvm1/celltypist/2ac906a5-9725-4258-8e36-21a9f6c0302a.rds")
meta <- typist@meta.data

subcl <- meta$Majority_voting_CellTypist
subcl[meta$tissue != "blood"] <- NA
cellgrp <- meta$Majority_voting_CellTypist_high
cellgrp[meta$tissue != "blood"] <- NA

mat <- typist@assays$RNA$counts

library(cellGeometry)
library(tictoc)

mk <- cellMarkers(typist@assays$RNA$counts, subclass = subcl,
                  cellgroup = cellgrp,
                  remove_subclass = c("Helper T cells", "Cytotoxic T cells"),
                  dual_mean = T)

toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog()

mk$mk_time <- mk_time

mk <- gene2symbol(mk, ensDb_v110)

#saveRDS(mk, "typist_cellmarkers_cycling.rds")

####generating simulated data####

library(cellGeometry)
library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

set.seed(3)
sim_counts <- generate_samples(mk, 25)
sim_percent <- sim_counts / rowSums(sim_counts) * 100
sim_pseudo <- simulate_bulk(mk, sim_counts)

sim_sampled_dir <- simulate_bulk(mat,
                                 sim_counts,
                                 subcl,
                                 method = "dirichlet") #times = 30 is default


rownames(sim_sampled_dir) <- gene2symbol(rownames(sim_sampled_dir), ensDb_v110)

# save(sim_counts, sim_percent, sim_pseudo,
#      sim_sampled_dir,
#      file = "simulated_dirichlet.rdata")

sim_sampled_dir_all <- list()

sim_sampled_dir_all[["Times 30"]][["Rep 1"]] <- sim_sampled_dir

for(x in paste("Rep", 2:5)){
  sim_sampled_dir_all[["Times 30"]][[x]] <- simulate_bulk(mat,
                                                          sim_counts,
                                                          subcl,
                                                          method = "dirichlet")
}

for(x in paste("Rep", 2:5)){
  rownames(sim_sampled_dir_all[["Times 30"]][[x]]) <- gene2symbol(rownames(sim_sampled_dir_all[["Times 30"]][[x]]), ensDb_v110)
}

# saveRDS(sim_sampled_dir_all, "simulated_dirichlet_all.rds")


####CellGeometry####

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
#   nsubclass  expfilter  use_filter  mean.RMSE
# 50        0.5       FALSE     0.6053

mk_update <- updateMarkers(mk, nsubclass = 50, expfilter = 0.5) 

geo_output <- function(times, rep){
  
  samples_sim <- sim_sampled_dir_all[[paste("Times", times)]][[paste("Rep", rep)]]
  tic()
  out <- deconvolute(mk_update, samples_sim,
                     arith_mean = T,
                     count_space = T, convert_bulk = F, 
                     use_filter = FALSE,
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
  cellgeo[["Times 30"]][[x]] <- geo_output(times = 30,
                                           rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(cellgeo, file = "cellgeo_dirichlet_output.rds")

####MuSiC2####

library(MuSiC)
library(SingleCellExperiment)

rownames(mat) <- gene2symbol(rownames(mat), ensDb_v110)

sce <- SingleCellExperiment(list(counts = mat[ , is.na(subcl) == FALSE]),
                            colData = data.frame("subclass" = subcl[is.na(subcl) == FALSE],
                                                 "donor_id" = meta$donor_id[is.na(subcl) == FALSE]))

music2_out <- function(times, rep){
  samples_sim <- sim_sampled_dir_all[[paste("Times", times)]][[paste("Rep", rep)]]
  
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

music2 <- list()

music2[["Times 30"]] <- list()

for(x in paste("Rep", 1:5)){
  music2[["Times 30"]][[x]] <- music2_out(times = 30,
                                          rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(music2, "music2_output.rds")

####DWLS####

# restart R/workstation 

typist <- readRDS("/media/lvm1/celltypist/2ac906a5-9725-4258-8e36-21a9f6c0302a.rds")

library(DWLS)
library(cellGeometry)
library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

#subset seurat 

sub <- subset(typist, subset = tissue == "blood")

meta <- sub@meta.data

subcl <- meta$Majority_voting_CellTypist
cellgrp <- meta$Majority_voting_CellTypist_high

dataSC  = sub@assays$RNA$counts
rownames(dataSC) <- gene2symbol(rownames(dataSC), ensDb_v110)
dataSC <- as.matrix(dataSC)

# problems with cluster names that DWLS accept and if there are 1 cell for subclass

subcl_tidy <- as.vector(subcl)
subcl_tidy[subcl_tidy == "Non-classical monocytes"] <- "Nonclassical monocytes"
subcl_tidy[subcl_tidy == "gamma-delta T cells"] <- "gamma_delta T cells"
subcl_tidy <- gsub("-", "neg", subcl_tidy)
subcl_tidy <- gsub("\\/", " ", subcl_tidy)
subcl_tidy <- gsub("\\+", "plus", subcl_tidy)
subcl_tidy <- gsub(" ", "_", subcl_tidy)

subcl_tidy <- factor(subcl_tidy)

remove <- names(table(subcl_tidy)[which(table(subcl_tidy) == 1)])
#"DC1"                    "Early erythroid"        "Type 17 helper T cells"

keep <- setdiff(unique(subcl_tidy), remove)

filter <- subcl_tidy %in% keep

subcl_tidyv2 <- subcl_tidy[filter]
subcl_tidyv2 <- factor(subcl_tidyv2)

dataSCv2 <- dataSC[ , filter]

meta <- meta[filter, ]
library(tictoc)
tic()
Signature <- buildSignatureMatrixMAST(scdata = dataSCv2,
                                      id = subcl_tidyv2,
                                      path = "DWLS_results")
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

# save(Signature, time, file = "DWLS_signature.rdata")

prop_arrangev2 <- function(df){
  colnames(df) <- gsub("_", " ", colnames(df))
  colnames(df) <- gsub("neg", "-", colnames(df))
  colnames(df)[colnames(df) == "Tem Effector helper T cells"] <- "Tem/Effector helper T cells"
  colnames(df)[colnames(df) == "Nonclassical monocytes"] <- "Non-classical monocytes"
  colnames(df)[colnames(df) == "CD16plus NK cells"] <- "CD16+ NK cells" 
  colnames(df)[colnames(df) == "Tcm Naive helper T cells"] <- "Tcm/Naive helper T cells" 
  colnames(df)[colnames(df) == "Tcm Naive cytotoxic T cells"] <- "Tcm/Naive cytotoxic T cells" 
  colnames(df)[colnames(df) == "Tem Effector cytotoxic T cells"] <- "Tem/Effector cytotoxic T cells" 
  colnames(df)[colnames(df) == "HSC MPP"] <- "HSC/MPP" 
  
  common <- intersect(colnames(sim_percent), colnames(df))
  df <- df[ , colnames(df) %in% common]
  df[ , match(colnames(sim_percent), colnames(df))]
}

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
    allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  }
  toc(log = TRUE) 
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  colnames(allCounts_DWLS) <- colnames(samples_sim)
  
  propdf <- prop_arrangev2(t(allCounts_DWLS * 100))
  
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

# saveRDS(DWLS, "DWLS_dirichlet_output.rds")

####LinDeconSeq####

# restart R/workstation 

typist <- readRDS("/media/lvm1/celltypist/2ac906a5-9725-4258-8e36-21a9f6c0302a.rds")

library(LinDeconSeq)
library(cellGeometry)
library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

#subset seurat 

sub <- subset(typist, subset = tissue == "blood")

meta <- sub@meta.data

subcl <- meta$Majority_voting_CellTypist
cellgrp <- meta$Majority_voting_CellTypist_high

dataSC  = sub@assays$RNA$counts
rownames(dataSC) <- gene2symbol(rownames(dataSC), ensDb_v110)
dataSC <- as.matrix(dataSC)

# again need to fix subclass names so it is acceptable
# remove subclass where 1 cell 

meta$cell_name <- rownames(meta)

remove <- names(table(subcl)[which(table(subcl) == 1)])
#"DC1"                    "Early erythroid"        "Type 17 helper T cells"

keep <- setdiff(unique(subcl), remove)

filter <- subcl %in% keep

table(filter)
#3 false

subclv2 <- subcl[filter]
subclv2 <- factor(subclv2)

dataSCv2 <- dataSC[ , filter]
rm(dataSC)

meta <- meta[filter, ]

meta$Majority_voting_CellTypist <- as.vector(meta$Majority_voting_CellTypist)
meta$Majority_voting_CellTypist[meta$Majority_voting_CellTypist == "Non-classical monocytes"] <- "Nonclassical monocytes"
meta$Majority_voting_CellTypist[meta$Majority_voting_CellTypist == "gamma-delta T cells"] <- "gamma_delta T cells"
meta$Majority_voting_CellTypist <- gsub("-", "neg", meta$Majority_voting_CellTypist)
meta$Majority_voting_CellTypist <- gsub("\\/", " ", meta$Majority_voting_CellTypist)
meta$Majority_voting_CellTypist <- gsub("\\+", "plus", meta$Majority_voting_CellTypist)
meta$Majority_voting_CellTypist <- gsub(" ", "_", meta$Majority_voting_CellTypist)

library(tidyr)
phes <- as.data.frame(meta %>% pivot_wider(Majority_voting_CellTypist,
                                           names_from = cell_name,
                                           values_from = cell_name))
rownames(phes) <- phes$Majority_voting_CellTypist
phes <- phes[ , -1]
phes <- as.matrix(phes)
phes <- ifelse(is.na(phes), 0, 1)

library(tictoc)
tic()
markerRes <- findMarkers(dataSCv2,
                         phes = phes)
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
markerRes$Time <- time

# saveRDS(markerRes, "LinDeconSeq_signature.rds")

prop_arrangev2 <- function(df){
  colnames(df) <- gsub("_", " ", colnames(df))
  colnames(df) <- gsub("neg", "-", colnames(df))
  colnames(df)[colnames(df) == "Tem Effector helper T cells"] <- "Tem/Effector helper T cells"
  colnames(df)[colnames(df) == "Nonclassical monocytes"] <- "Non-classical monocytes"
  colnames(df)[colnames(df) == "CD16plus NK cells"] <- "CD16+ NK cells" 
  colnames(df)[colnames(df) == "Tcm Naive helper T cells"] <- "Tcm/Naive helper T cells" 
  colnames(df)[colnames(df) == "Tcm Naive cytotoxic T cells"] <- "Tcm/Naive cytotoxic T cells" 
  colnames(df)[colnames(df) == "Tem Effector cytotoxic T cells"] <- "Tem/Effector cytotoxic T cells" 
  colnames(df)[colnames(df) == "HSC MPP"] <- "HSC/MPP" 
  
  common <- intersect(colnames(sim_percent), colnames(df))
  df <- df[ , colnames(df) %in% common]
  df[ , match(colnames(sim_percent), colnames(df))]
}

linDecon_out <- function(times, rep){
  samples_sim <- sim_sampled_dir_all[[paste("Times", times)]][[paste("Rep", rep)]]
  
  tic()
  fractions <- deconSeq(samples_sim , markerRes$sigMatrix$sig.mat, verbose = FALSE)
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  propdf <- prop_arrangev2(fractions * 100)
  
  metric_percent <- metric_set(sim_percent, propdf)
  
  list(output = fractions,
       metrics_percent = metric_percent,
       time = time)
  
}

Lin <- list()

Lin[["Times 30"]] <- list()

for(x in paste("Rep", 1:5)){
  Lin[["Times 30"]][[x]] <- linDecon_out(times = 30,
                                         rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(Lin, "LinDeconSeq_output.rds")






