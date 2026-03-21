####Background####

# Assessing NMF deconvolution in AMP simulated data

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet")

load("/media/gcpeac/Myles/Deconvolution/AMP_scRNAseq_data_and_annotations.RData")

celseq_counts <- as.matrix(celseq_counts)

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

sim_sampled_dir_all <- readRDS("simulated_dirichlet_all.rds")
load("simulated_dirichlet.rdata")

library(cellGeometry)
mk <- readRDS("cellmarkers_new.rds")

####NMF####

library(NMF)

mkx <- 2^mk$genemeans_ar -1

targ <- sim_sampled_dir_all$`Times 30`$`Rep 1`[rownames(mkx), ]
mode(targ) <- "numeric"

#all(rownames(targ) == rownames(mkx))
#TRUE

ok <- rowSums(sim_sampled_dir_all$`Times 30`$`Rep 1`[rownames(mkx), ]) > 0
table(ok)
# FALSE  TRUE 
# 569 24508 

# test first then function

# Initial factorisation 
init <- nmfModel(rank = ncol(mkx), target = targ[ok, ], W = mkx[ok, ], H = 1)
str(init)

library(tictoc)
tic()
res <- nmf(targ[ok, ], ncol(mkx), seed = init) #by default it uses brunet 
toc()

str(res)
#res@fit@H

pred_counts <- res@fit@H
pred_percent <- apply(pred_counts, 2, function(x){x/sum(x)})
pred_percent <- pred_percent * 100
pred_percent <- t(pred_percent)

#plot_set(sim_percent, pred_percent)

rm(targ)

nmf_out <- function(df, times){
  targ <- df[rownames(mkx), ]
  mode(targ) <- "numeric"
  
  ok <- rowSums(df[rownames(mkx), ]) > 0
  tic()
  init <- nmfModel(rank = ncol(mkx), target = targ[ok, ], W = mkx[ok, ], H = 1)
  res <- nmf(targ[ok, ], ncol(mkx), seed = init) 
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  pred_counts <- res@fit@H
  pred_percent <- apply(pred_counts, 2, function(x){x/sum(x)})
  pred_percent <- pred_percent * 100
  pred_percent <- t(pred_percent)
  pred_percent <- pred_percent[ , colnames(sim_percent)]
  
  pred_counts <- t(pred_counts)
  pred_counts <- pred_counts[ , colnames(sim_counts)]
  
  metrics_counts <- metric_set(sim_counts, pred_counts/times)
  metrics_percent <- metric_set(sim_percent, pred_percent)
  
  list(output = res,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = time)
}

AMP_NMF <- list()

for(x in paste("Rep", 1:5)){
  AMP_NMF[["Times 30"]][[x]] <- nmf_out(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                        times = 30)
}

#saveRDS(AMP_NMF, "NMF_dirichlet_output.rds")

#identical(AMP_NMF$`Times 30`$`Rep 1`$output@fit@H, res@fit@H)
#TRUE

#####other NMF methods#####

meths <- c("KL", "lee", "snmf/r")

mkx <- 2^mk$genemeans_ar -1

targ <- sim_sampled_dir_all$`Times 30`$`Rep 1`[rownames(mkx), ]
mode(targ) <- "numeric"

#all(rownames(targ) == rownames(mkx))
#TRUE

ok <- rowSums(sim_sampled_dir_all$`Times 30`$`Rep 1`[rownames(mkx), ]) > 0
table(ok)
# FALSE  TRUE 
# 569 24508 

init <- nmfModel(rank = ncol(mkx), target = targ[ok, ], W = mkx[ok, ], H = 1)

library(pbmclapply)

res <- lapply(meths, function(i) {
  out <- nmf(targ[ok, ], ncol(mkx), seed = init, method = i)
  out
})

names(res) <- meths

res$`snmf/r`@fit@H
#NMF::snmf - Too many restarts due to too big 'beta' value 
#Computation stopped after the 9th restart

out <- nmf(targ[ok, ], ncol(mkx), seed = init, method = "snmf/l")
#nsNMF error
#pe-nmf error
#snmf/l error

rm(targ)

nmf_out_alt <- function(df, times, method){
  targ <- df[rownames(mkx), ]
  mode(targ) <- "numeric"
  
  ok <- rowSums(df[rownames(mkx), ]) > 0
  tic()
  init <- nmfModel(rank = ncol(mkx), target = targ[ok, ], W = mkx[ok, ], H = 1)
  res <- nmf(targ[ok, ], ncol(mkx), seed = init, method = method) 
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  pred_counts <- res@fit@H
  pred_percent <- apply(pred_counts, 2, function(x){x/sum(x)})
  pred_percent <- pred_percent * 100
  pred_percent <- t(pred_percent)
  pred_percent <- pred_percent[ , colnames(sim_percent)]
  
  pred_counts <- t(pred_counts)
  pred_counts <- pred_counts[ , colnames(sim_counts)]
  
  metrics_counts <- metric_set(sim_counts, pred_counts/times)
  metrics_percent <- metric_set(sim_percent, pred_percent)
  
  list(output = res,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = time)
}

AMP_NMF_alt <- list()

for(i in meths){
  for(x in paste("Rep", 1:5)){
    AMP_NMF_alt[[i]][[x]] <- nmf_out_alt(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                              times = 30,
                                              method = i)
  }
}

#nmf_snmf failed to run

#saveRDS(AMP_NMF_alt, "NMF_dirichlet_output_alt.rds")

#restart R

#RcppML
library(RcppML)
library(cellGeometry)
library(tictoc)

out <- project(targ[ok, ], w = mkx[ok, ])

RcppML_out <- function(df, times){
  targ <- df[rownames(mkx), ]
  mode(targ) <- "numeric"
  
  ok <- rowSums(df[rownames(mkx), ]) > 0
  tic()
  out <- project(targ[ok, ], w = mkx[ok, ])
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  pred_counts <- out
  rownames(pred_counts) <- colnames(mkx)
  pred_percent <- apply(pred_counts, 2, function(x){x/sum(x)})
  pred_percent <- pred_percent * 100
  pred_percent <- t(pred_percent)
  pred_percent <- pred_percent[ , colnames(sim_percent)]
  
  pred_counts <- t(pred_counts)
  pred_counts <- pred_counts[ , colnames(sim_counts)]
  
  metrics_counts <- metric_set(sim_counts, pred_counts/times)
  metrics_percent <- metric_set(sim_percent, pred_percent)
  
  list(output = out,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = time)
}

AMP_RcppML <- list()


for(x in paste("Rep", 1:5)){
  AMP_RcppML[["Times 30"]][[x]] <- RcppML_out(df = sim_sampled_dir_all[["Times 30"]][[x]],
                                       times = 30)
}


#saveRDS(AMP_RcppML, "RcppML_dirichlet_output.rds")