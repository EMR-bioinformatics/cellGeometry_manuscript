####Background####

# Accessing NMF deconvolution for brain dataset

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain")

load("brain_simulated_merged.rdata")

mk <- readRDS("brain_cellmarkers_merged.rds")

####NMF####

library(NMF)

mkx <- 2^mk$genemeans_ar -1

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
  pred_percent <- pred_percent[ , colnames(sim_percent_merge)]
  
  pred_counts <- t(pred_counts)
  pred_counts <- pred_counts[ , colnames(sim_counts_merge)]
  
  metrics_counts <- metric_set(sim_counts_merge, pred_counts/times)
  metrics_percent <- metric_set(sim_percent_merge, pred_percent)
  
  list(output = res,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = time)
}

NMF <- list()

library(tictoc)

for(x in paste("Rep", 1:5)){
  NMF[["Times 1"]][[x]] <- nmf_out(df = sim_sampled_merge[["Times 1"]][[x]],
                                   times = 3)
}

#saveRDS(NMF, "NMF_dirichlet_output.rds")

#####other NMF methods#####

meths <- c("KL", "lee")

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
  pred_percent <- pred_percent[ , colnames(sim_percent_merge)]
  
  pred_counts <- t(pred_counts)
  pred_counts <- pred_counts[ , colnames(sim_counts_merge)]
  
  metrics_counts <- metric_set(sim_counts_merge, pred_counts/times)
  metrics_percent <- metric_set(sim_percent_merge, pred_percent)
  
  list(output = res,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = time)
}

NMF_alt <- list()

for(i in meths){
  for(x in paste("Rep", 1:5)){
    NMF_alt[[i]][[x]] <- nmf_out_alt(df = sim_sampled_merge[["Times 1"]][[x]],
                                     times = 1,
                                     method = i)
  }
}

#saveRDS(NMF_alt, "NMF_dirichlet_output_alt.rds")
# Warning messages:
#   1: In validityMethod(object) :
#   Dimensions of W and H look strange [ncol(W)= 116 > ncol(H)= 30 ]

#Try the other NMF methods
targ <- sim_sampled_merge$`Times 1`$`Rep 1`[rownames(mkx), ]
mode(targ) <- "numeric"

#all(rownames(targ) == rownames(mkx))
#TRUE

ok <- rowSums(sim_sampled_merge$`Times 1`$`Rep 1`[rownames(mkx), ]) > 0
table(ok)
# FALSE  TRUE 
# 1610 60149 

init <- nmfModel(rank = ncol(mkx), target = targ[ok, ], W = mkx[ok, ], H = 1)

out <- nmf(targ[ok, ], ncol(mkx), seed = init, method = "snmf/l")
# NMF::snmf - Too many restarts due to too big 'beta' value [Computation stopped after the 9th restart]
# nsNMF error
# pe-nmf error
# snmf/l error

# restart R

library(RcppML)
library(cellGeometry)
library(tictoc)

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
  pred_percent <- pred_percent[ , colnames(sim_percent_merge)]
  
  pred_counts <- t(pred_counts)
  pred_counts <- pred_counts[ , colnames(sim_counts_merge)]
  
  metrics_counts <- metric_set(sim_counts_merge, pred_counts/times)
  metrics_percent <- metric_set(sim_percent_merge, pred_percent)
  
  list(output = out,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = time)
}

RcppML <- list()

for(x in paste("Rep", 1:5)){
  RcppML[["Times 1"]][[x]] <- RcppML_out(df = sim_sampled_merge[["Times 1"]][[x]],
                                         times = 3)
}


#saveRDS(RcppML, "RcppML_dirichlet_output.rds")



