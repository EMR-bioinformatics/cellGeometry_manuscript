####Background####

# Benchmarking of glmnet with various datasets using cellGeometry gene signature 

####AMP####

library(cellGeometry)
mk <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/cellmarkers_new.rds")
mkx <- 2^mk$genemeans_ar -1

sim_sampled_dir_all <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/simulated_dirichlet_all.rds")
load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/simulated_dirichlet.rdata")

library(glmnet)

ok <- mk$geneset

genes <- intersect(rownames(mkx), rownames(exp))

exp <- sim_sampled_dir_all$`Times 30`$`Rep 1`

library(mcprogress)

res <- pmclapply(seq_len(ncol(exp)), function(i) {
  fit <- cv.glmnet(x = mkx[genes, ], y = exp[genes, i], lower.limits = 0,
                   alpha = 1)
  cf <- coef(fit)
  cf[-1, 1]
}, mc.cores = 8)

res <- do.call(cbind, res)
res <- t(res)
colnames(res) <- colnames(mkx)
res_percent <- apply(res, 1, function(x){x/sum(x)})
res_percent <- t(res_percent)

library(tictoc)

Rsq <- function(obs, pred) {
  rss <- sum((pred - obs)^2)
  tss <- sum((obs - mean(obs))^2)
  1 - rss/tss
}

glmnet_out <- function(df, times, alpha){
  
  tic()
  temp <- pmclapply(seq_len(ncol(df)), function(i){
    fit <- cv.glmnet(x = mkx[genes, ], y = df[genes, i], lower.limits = 0,
                     alpha = alpha)
    cf <- coef(fit)
    cf[-1, 1]
  }, mc.cores = 8) 
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
       rsq_counts = Rsq(sim_counts, counts/times),
       rsq_percent = Rsq(sim_percent, percent),
       time = time)
  
}

amp_glmnet <- list()

for(i in seq(0, 1, 0.1)){
  
  amp_glmnet[["Times 30"]][[paste0("alpha = ", i)]] <- list()
  
  for(x in paste("Rep", 1:5)){
    amp_glmnet[["Times 30"]][[paste0("alpha = ", i)]][[x]] <- glmnet_out(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                                                         times = 30,
                                                                         alpha = i)
  }
  
}

# saveRDS(amp_glmnet, file = "glmnet_amp_full.rds")

glmnet_out_sub <- function(df, times, alpha){
  
  tic()
  temp <- lapply(seq_len(ncol(df)), function(i){
    fit <- cv.glmnet(x = mkx[ok, ], y = df[ok, i], lower.limits = 0,
                     alpha = alpha)
   
    cf <- coef(fit)
    cf[-1, 1]
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
       rsq_counts = Rsq(sim_counts, counts/times),
       rsq_percent = Rsq(sim_percent, percent),
       time = time)
  
}

amp_glmnet_sub <- list()

for(i in seq(0, 1, 0.1)){
  
  amp_glmnet_sub[["Times 30"]][[paste0("alpha = ", i)]] <- list()
  
  for(x in paste("Rep", 1:5)){
    amp_glmnet_sub[["Times 30"]][[paste0("alpha = ", i)]][[x]] <- glmnet_out_sub(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                                                                 times = 30,
                                                                                 alpha = i)
  }
  
}

# saveRDS(amp_glmnet_sub, file = "glmnet_amp.rds") 

####typist####

rm(list = setdiff(ls(), c("glmnet_out", "glmnet_out_sub", "Rsq")))

path <- "/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Cell_typist/workstation2/Dirichlet"

mk <- readRDS(paste0(path, "/typist_cellmarkers_cycling.rds"))
sim_sampled_dir_all <- readRDS(paste0(path, "/simulated_dirichlet_all.rds"))
load(paste0(path, "/simulated_dirichlet.rdata"))

mkx <- 2^mk$genemeans_ar -1

ok <- mk$geneset

genes <- intersect(rownames(mkx), rownames(sim_sampled_dir_all$`Times 30`$`Rep 1`))

typist_glmnet <- list()

for(i in seq(0, 1, 0.1)){
  
  typist_glmnet[["Times 30"]][[paste0("alpha = ", i)]] <- list()
  
  for(x in paste("Rep", 1:5)){
    typist_glmnet[["Times 30"]][[paste0("alpha = ", i)]][[x]] <- glmnet_out(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                                                         times = 30,
                                                                         alpha = i)
  }
  
}

#saveRDS(typist_glmnet, file = "glmnet_typist_full.rds")

typist_glmnet_sub <- list()

for(i in seq(0, 1, 0.1)){
  
  typist_glmnet_sub[["Times 30"]][[paste0("alpha = ", i)]] <- list()
  
  for(x in paste("Rep", 1:5)){
    typist_glmnet_sub[["Times 30"]][[paste0("alpha = ", i)]][[x]] <- glmnet_out_sub(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                                                                 times = 30,
                                                                                 alpha = i)
  }
  
}

#saveRDS(typist_glmnet_sub, file = "glmnet_typist.rds") 

####tabula####

rm(list = setdiff(ls(), c("glmnet_out", "glmnet_out_sub", "Rsq")))

path <- "/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2"

load(paste0(path, "/tabula_simulated_workstation2.rdata"))

mk <- readRDS(paste0(path, "/tabula_markers_dualmeans_rerun.rds"))

mkx <- 2^mk$genemeans_ar -1

ok <- mk$geneset

genes <- intersect(rownames(mkx), rownames(sim_sampled$`Times 3`$`Rep 1`))

tabula_glmnet <- list()

for(i in seq(0, 1, 0.1)){
  
  tabula_glmnet[["Times 3"]][[paste0("alpha = ", i)]] <- list()
  
  for(x in paste("Rep", 1:5)){
    tabula_glmnet[["Times 3"]][[paste0("alpha = ", i)]][[x]] <- glmnet_out(df = sim_sampled[["Times 3"]][[x]],
                                                                            times = 3,
                                                                            alpha = i)
  }
  
}

#saveRDS(tabula_glmnet, file = "glmnet_tabula_full.rds")
#getting long 
#will focus on alpha 0 and 1 for bigger datasets

tabula_glmnet_sub <- list()

for(i in c(0, 1)){
  
  tabula_glmnet_sub[["Times 3"]][[paste0("alpha = ", i)]] <- list()
  
  for(x in paste("Rep", 1:5)){
    tabula_glmnet_sub[["Times 3"]][[paste0("alpha = ", i)]][[x]] <- glmnet_out_sub(df = sim_sampled[["Times 3"]][[x]],
                                                                                    times = 3,
                                                                                    alpha = i)
  }
  
}

# only do 0 and 1

#saveRDS(tabula_glmnet_sub, file = "glmnet_tabula.rds")

####brain####

rm(list = setdiff(ls(), "Rsq"))

path <- "/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain"

load(paste0(path, "/brain_simulated_merged.rdata"))

mk <- readRDS(paste0(path, "/brain_cellmarkers_merged.rds"))

mkx <- 2^mk$genemeans_ar -1

ok <- mk$geneset

genes <- intersect(rownames(mkx), rownames(sim_sampled_merge$`Times 1`$`Rep 1`))

glmnet_out <- function(df, times, alpha){
  
  tic()
  temp <- pmclapply(seq_len(ncol(df)), function(i){
    fit <- cv.glmnet(x = mkx[genes, ], y = df[genes, i], lower.limits = 0,
                     alpha = alpha)
    cf <- coef(fit)
    cf[-1, 1]
  }, mc.cores = 8) 
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
       rsq_counts = Rsq(sim_counts_merge, counts/times),
       rsq_percent = Rsq(sim_percent_merge, percent),
       time = time)
  
}

brain_glmnet <- list()

for(i in c(0, 1)){
  
  brain_glmnet[["Times 1"]][[paste0("alpha = ", i)]] <- list()
  
  for(x in paste("Rep", 1:5)){
    brain_glmnet[["Times 1"]][[paste0("alpha = ", i)]][[x]] <- glmnet_out(df = sim_sampled_merge[["Times 1"]][[x]],
                                                                           times = 1,
                                                                           alpha = i)
  }
  
}

#saveRDS(brain_glmnet, file = "glmnet_brain_full.rds")

glmnet_out_sub <- function(df, times, alpha){
  
  tic()
  temp <- lapply(seq_len(ncol(df)), function(i){
    fit <- cv.glmnet(x = mkx[ok, ], y = df[ok, i], lower.limits = 0,
                     alpha = alpha)
    
    cf <- coef(fit)
    cf[-1, 1]
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
       rsq_counts = Rsq(sim_counts_merge, counts/times),
       rsq_percent = Rsq(sim_percent_merge, percent),
       time = time)
  
}

brain_glmnet_sub <- list()

for(i in c(0, 1)){
  
  brain_glmnet_sub[["Times 1"]][[paste0("alpha = ", i)]] <- list()
  
  for(x in paste("Rep", 1:5)){
    brain_glmnet_sub[["Times 1"]][[paste0("alpha = ", i)]][[x]] <- glmnet_out_sub(df = sim_sampled_merge[["Times 1"]][[x]],
                                                                                   times = 3,
                                                                                   alpha = i)
  }
  
}

#saveRDS(brain_glmnet_sub, file = "glmnet_brain.rds")

