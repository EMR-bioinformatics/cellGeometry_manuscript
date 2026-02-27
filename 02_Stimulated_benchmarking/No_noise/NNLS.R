####Background####

# Benchmarking of NNLS with various datasets using cellGeometry gene signature

####AMP####

library(cellGeometry)
mk <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/cellmarkers_new.rds")
mkx <- 2^mk$genemeans_ar -1

sim_sampled_dir_all <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/simulated_dirichlet_all.rds")
load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/simulated_dirichlet.rdata")

library(nnls)

ok <- mk$geneset

exp <- sim_sampled_dir_all$`Times 30`$`Rep 1`

res <- lapply(seq_len(ncol(exp)), function(i) {
  fit <- nnls(A = mkx[ok, ], b = exp[ok, i])
  fit$x
})

res <- do.call(cbind, res)
res <- t(res)
colnames(res) <- colnames(mkx)
res_percent <- apply(res, 1, function(x){x/sum(x)})
res_percent <- t(res_percent)

library(tictoc)

nnls_out <- function(df, times){
  
  tic()
  temp <- lapply(seq_len(ncol(df)), function(i) {
    fit <- nnls(A = mkx[ok, ], b = df[ok, i])
    fit$x
  })
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  counts <- do.call(cbind, temp)
  counts <- t(counts)
  colnames(counts) <- colnames(mkx)
  
  percent <- apply(counts, 1, function(x){x/sum(x)})
  percent <- t(percent)
  
  metrics_counts <- metric_set(sim_counts, counts/times)
  metrics_percent <- metric_set(sim_percent, percent)
  
  list(output = counts,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = time)
}

AMP_nnls <- list()

for(x in paste("Rep", 1:5)){
  AMP_nnls[["Times 30"]][[x]] <- nnls_out(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                         times = 30)
}

#saveRDS(AMP_nnls, "nnls_AMP.rds")

####typist####

rm(list = setdiff(ls(), "nnls_out"))

path <- "/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Cell_typist/workstation2/Dirichlet"

mk <- readRDS(paste0(path, "/typist_cellmarkers_cycling.rds"))
sim_sampled_dir_all <- readRDS(paste0(path, "/simulated_dirichlet_all.rds"))
load(paste0(path, "/simulated_dirichlet.rdata"))

mkx <- 2^mk$genemeans_ar -1

ok <- mk$geneset

typist_nnls <- list()

for(x in paste("Rep", 1:5)){
  typist_nnls[["Times 30"]][[x]] <- nnls_out(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                          times = 30)
}

#saveRDS(typist_nnls, "nnls_typist.rds")

####tabula####

rm(list = setdiff(ls(), "nnls_out"))

path <- "/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2"

load(paste0(path, "/tabula_simulated_workstation2.rdata"))

mk <- readRDS(paste0(path, "/tabula_markers_dualmeans_rerun.rds"))

mkx <- 2^mk$genemeans_ar -1

ok <- mk$geneset

tabula_nnls <- list()

for(x in paste("Rep", 1:5)){
  tabula_nnls[["Times 3"]][[x]] <- nnls_out(df = sim_sampled[["Times 3"]][[x]],
                                             times = 3)
}

#saveRDS(tabula_nnls, "nnls_tabula.rds")

####brain####

rm(list = ls())

path <- "/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain"

load(paste0(path, "/brain_simulated_merged.rdata"))

mk <- readRDS(paste0(path, "/brain_cellmarkers_merged.rds"))

mkx <- 2^mk$genemeans_ar -1

ok <- mk$geneset

nnls_out <- function(df, times){
  
  tic()
  temp <- lapply(seq_len(ncol(df)), function(i) {
    fit <- nnls(A = mkx[ok, ], b = df[ok, i])
    fit$x
  })
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  counts <- do.call(cbind, temp)
  counts <- t(counts)
  colnames(counts) <- colnames(mkx)
  
  percent <- apply(counts, 1, function(x){x/sum(x)})
  percent <- t(percent)
  
  metrics_counts <- metric_set(sim_counts_merge, counts/times)
  metrics_percent <- metric_set(sim_percent_merge, percent)
  
  list(output = counts,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = time)
}

brain_nnls <- list()

for(x in paste("Rep", 1:5)){
  brain_nnls[["Times 1"]][[x]] <- nnls_out(df = sim_sampled_merge[["Times 1"]][[x]],
                                            times = 1)
}

#saveRDS(brain_nnls, "nnls_brain.rds")





