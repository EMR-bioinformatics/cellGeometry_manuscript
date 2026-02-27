####Background####

# Different proportion of cell subclasses set to zero 
# workstation 1

load("/media/gcpeac/Myles/Deconvolution/AMP_scRNAseq_data_and_annotations.RData")

celseq_counts <- as.matrix(celseq_counts)

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

####generate samples + bulk####

library(cellGeometry)

mk <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/cellmarkers_new.rds")

sim_counts <- generate_samples(mk, 25, zero_fraction = 0.1)
sim_percent <- sim_counts / rowSums(sim_counts) * 100

sim_counts20 <- generate_samples(mk, 25, zero_fraction = 0.2)
sim_percent20 <- sim_counts20 / rowSums(sim_counts20) * 100

sim_sampled_dir_all <- list()

sim_sampled_dir_all[["zero_10"]] <- list()
sim_sampled_dir_all[["zero_10"]][["Times 30"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampled_dir_all[["zero_10"]][["Times 30"]][[x]] <- simulate_bulk(celseq_counts,
                                                                       sim_counts,
                                                                       metadata$subclass,
                                                                       method = "dirichlet")
}

sim_sampled_dir_all[["zero_20"]] <- list()
sim_sampled_dir_all[["zero_20"]][["Times 30"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampled_dir_all[["zero_20"]][["Times 30"]][[x]] <- simulate_bulk(celseq_counts,
                                                                       sim_counts20,
                                                                       metadata$subclass,
                                                                       method = "dirichlet")
}

# save(sim_counts, sim_counts20,
#      sim_percent, sim_percent20,
#      sim_sampled_dir_all,
#      file = "simulated_dirichlet_zerofraction.rdata")

#generate other fractions
load("simulated_dirichlet_zerofraction.rdata")

sim_counts_none <- generate_samples(mk, 25, zero_fraction = 0)
sim_percent_none <- sim_counts_none / rowSums(sim_counts_none) * 100

sim_counts5 <- generate_samples(mk, 25, zero_fraction = 0.05)
sim_percent5 <- sim_counts5 / rowSums(sim_counts5) * 100

sim_counts50 <- generate_samples(mk, 25, zero_fraction = 0.5)
sim_percent50 <- sim_counts50 / rowSums(sim_counts50) * 100

sim_sampled_dir_all[["zero_0"]] <- list()
sim_sampled_dir_all[["zero_0"]][["Times 30"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampled_dir_all[["zero_0"]][["Times 30"]][[x]] <- simulate_bulk(celseq_counts,
                                                                       sim_counts_none,
                                                                       metadata$subclass,
                                                                       method = "dirichlet")
}

sim_sampled_dir_all[["zero_5"]] <- list()
sim_sampled_dir_all[["zero_5"]][["Times 30"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampled_dir_all[["zero_5"]][["Times 30"]][[x]] <- simulate_bulk(celseq_counts,
                                                                      sim_counts5,
                                                                      metadata$subclass,
                                                                      method = "dirichlet")
}

sim_sampled_dir_all[["zero_50"]] <- list()
sim_sampled_dir_all[["zero_50"]][["Times 30"]] <- list()

for(x in paste("Rep", 1:5)){
  sim_sampled_dir_all[["zero_50"]][["Times 30"]][[x]] <- simulate_bulk(celseq_counts,
                                                                      sim_counts50,
                                                                      metadata$subclass,
                                                                      method = "dirichlet")
}

# save(sim_counts, sim_counts20,
#      sim_percent, sim_percent20,
#      sim_sampled_dir_all,
#      sim_counts_none, sim_percent_none,
#      sim_counts5, sim_percent5,
#      sim_counts50, sim_percent50,
#      file = "simulated_dirichlet_zerofraction.rdata")

####cellGeometry####

mk_update <- updateMarkers(mk, nsubclass = 500, expfilter = 0.25)
#same parameter as AMP no zero fraction 

library(tictoc)
geo_output <- function(df, times, counts, percent){
  
  tic()
  out <- deconvolute(mk_update, df,
                     arith_mean = T,
                     count_space = T, convert_bulk = F, 
                     use_filter = FALSE,
                     weight_method = "equal")

  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  metrics_counts <- metric_set(counts, out$subclass$output/times)
  metrics_percent <- metric_set(percent, out$subclass$percent)
  
  list(output = out,
       metrics_counts = metrics_counts,
       metrics_percent = metrics_percent,
       time = mk_time)
  
}

cellgeo <- list()

for(x in paste("Rep", 1:5)){
  cellgeo[["zero_10"]][["Times 30"]][[x]] <- geo_output(sim_sampled_dir_all$zero_10$`Times 30`[[x]],
                                                        times = 30, sim_counts, sim_percent)
}

for(x in paste("Rep", 1:5)){
  cellgeo[["zero_20"]][["Times 30"]][[x]] <- geo_output(sim_sampled_dir_all$zero_20$`Times 30`[[x]],
                                                        times = 30, sim_counts20, sim_percent20)
}

#saveRDS(cellgeo, file = "cellgeo_zerofraction.rds")

cellgeo <- readRDS("cellgeo_zerofraction.rds")

for(x in paste("Rep", 1:5)){
  cellgeo[["zero_0"]][["Times 30"]][[x]] <- geo_output(sim_sampled_dir_all$zero_0$`Times 30`[[x]],
                                                        times = 30, sim_counts_none, sim_percent_none)
}

for(x in paste("Rep", 1:5)){
  cellgeo[["zero_5"]][["Times 30"]][[x]] <- geo_output(sim_sampled_dir_all$zero_5$`Times 30`[[x]],
                                                       times = 30, sim_counts5, sim_percent5)
}

for(x in paste("Rep", 1:5)){
  cellgeo[["zero_50"]][["Times 30"]][[x]] <- geo_output(sim_sampled_dir_all$zero_50$`Times 30`[[x]],
                                                       times = 30, sim_counts50, sim_percent50)
}

#saveRDS(cellgeo, file = "cellgeo_zerofraction.rds")

####MuSiC####

library(SingleCellExperiment)
library(MuSiC, lib.loc =  "/home/rachel/R/x86_64-pc-linux-gnu-library/4.2")

sce <- SingleCellExperiment(list(counts = celseq_counts[ ,is.na(metadata$subclass) == FALSE]),
                            colData = metadata[is.na(metadata$subclass) == FALSE, ])

library(tictoc)

music2_out <- function(df, percent){
  
  tic()
  out <- music_prop(bulk.mtx = df,
                    sc.sce = sce,
                    clusters = "subclass",
                    samples = "sample")
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  propdf <- out$Est.prop.weighted * 100
  propdf <- propdf[ , match(colnames(percent), colnames(propdf))]
  
  metrics_percent <- metric_set(percent, propdf)
  
  list(output = out,
       metrics_percent = metrics_percent,
       time = mk_time)
}

music_res <- list()

for(x in paste("Rep", 1:5)){
  music_res[["zero_10"]][["Times 30"]][[x]] <- music2_out(sim_sampled_dir_all$zero_10$`Times 30`[[x]],
                                                          sim_percent)
}

for(x in paste("Rep", 1:5)){
  music_res[["zero_20"]][["Times 30"]][[x]] <- music2_out(sim_sampled_dir_all$zero_20$`Times 30`[[x]],
                                                          sim_percent20)
}


#saveRDS(music_res, file = "music_zerofraction.rds")

music_res <- readRDS("music_zerofraction.rds")

for(x in paste("Rep", 1:5)){
  music_res[["zero_0"]][["Times 30"]][[x]] <- music2_out(sim_sampled_dir_all$zero_0$`Times 30`[[x]],
                                                          sim_percent_none)
}

for(x in paste("Rep", 1:5)){
  music_res[["zero_5"]][["Times 30"]][[x]] <- music2_out(sim_sampled_dir_all$zero_5$`Times 30`[[x]],
                                                         sim_percent5)
}

for(x in paste("Rep", 1:5)){
  music_res[["zero_50"]][["Times 30"]][[x]] <- music2_out(sim_sampled_dir_all$zero_50$`Times 30`[[x]],
                                                         sim_percent50)
}

#saveRDS(music_res, file = "music_zerofraction.rds")

####DWLS####

library(DWLS)
load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Exploring_other_tools/DWLS_test/Signature.rdata")
rm(time)

DWLS_out <- function(df, percent){
  allCounts_DWLS<-NULL
  
  tic() 
  for(j in 1:(dim(df)[2])){
    S <- Signature
    Bulk<-df[, j]
    names(Bulk)<- rownames(df)
    Genes <- intersect(rownames(S), names(Bulk))
    B <- Bulk[Genes]
    S <- S[Genes,]
    solDWLS<-try(solveDampenedWLS(S,B), silent = TRUE)
    if (inherits(solDWLS, 'try-error')) solDWLS <- rep(NA, ncol(S))
    #try is a new addition since had problems with times = 1, rep 4
    allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  }
  toc(log = TRUE) 
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  colnames(allCounts_DWLS) <- colnames(df)
  colnames(allCounts_DWLS) <- colnames(df)
  rownames(allCounts_DWLS) <- paste0("SC-", rownames(allCounts_DWLS))
  propdf <- t(allCounts_DWLS) * 100
  propdf <- propdf[ , match(colnames(percent),
                            colnames(propdf))] 
  metric_percent <- metric_set(percent, propdf)
  
  list(output = allCounts_DWLS,
       metrics_percent = metric_percent,
       time = time)
  
}

DWLS_res <- list()

for(x in paste("Rep", 1:5)){
  DWLS_res[["zero_10"]][["Times 30"]][[x]] <- DWLS_out(sim_sampled_dir_all$zero_10$`Times 30`[[x]],
                                                          sim_percent)
}

for(x in paste("Rep", 1:5)){
  DWLS_res[["zero_20"]][["Times 30"]][[x]] <- DWLS_out(sim_sampled_dir_all$zero_20$`Times 30`[[x]],
                                                          sim_percent20)
}


#saveRDS(DWLS_res, file = "DWLS_zerofraction.rds")

DWLS_res <- readRDS("DWLS_zerofraction.rds")

for(x in paste("Rep", 1:5)){
  DWLS_res[["zero_0"]][["Times 30"]][[x]] <- DWLS_out(sim_sampled_dir_all$zero_0$`Times 30`[[x]],
                                                       sim_percent_none)
}

for(x in paste("Rep", 1:5)){
  DWLS_res[["zero_5"]][["Times 30"]][[x]] <- DWLS_out(sim_sampled_dir_all$zero_5$`Times 30`[[x]],
                                                      sim_percent5)
}

for(x in paste("Rep", 1:5)){
  DWLS_res[["zero_50"]][["Times 30"]][[x]] <- DWLS_out(sim_sampled_dir_all$zero_50$`Times 30`[[x]],
                                                      sim_percent50)
}

#saveRDS(DWLS_res, file = "DWLS_zerofraction.rds")


####LinDeconSeq####

# restart R

markerRes <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/LinDeconSeq_signature.rds")

library(LinDeconSeq)
library(cellGeometry)
library(tictoc)

linDecon_out <- function(df, percent){
  
  tic()
  fractions <- deconSeq(df, markerRes$sigMatrix$sig.mat, verbose = FALSE)
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  propdf <- fractions * 100
  colnames(propdf) <- gsub("\\.", "-", colnames(propdf))
  propdf <- propdf[ , match(colnames(percent),
                            colnames(propdf))]
  
  metric_percent <- metric_set(percent, propdf)
  
  list(output = fractions,
       metrics_percent = metric_percent,
       time = time)
  
}

lin <- list()

for(x in paste("Rep", 1:5)){
  lin[["zero_10"]][["Times 30"]][[x]] <- linDecon_out(sim_sampled_dir_all$zero_10$`Times 30`[[x]],
                                                       sim_percent)
}

for(x in paste("Rep", 1:5)){
  lin[["zero_20"]][["Times 30"]][[x]] <- linDecon_out(sim_sampled_dir_all$zero_20$`Times 30`[[x]],
                                                       sim_percent20)
}


#saveRDS(lin, file = "lin_zerofraction.rds")

lin <- readRDS("lin_zerofraction.rds")

for(x in paste("Rep", 1:5)){
  lin[["zero_0"]][["Times 30"]][[x]] <- linDecon_out(sim_sampled_dir_all$zero_0$`Times 30`[[x]],
                                                      sim_percent_none)
}

for(x in paste("Rep", 1:5)){
  lin[["zero_5"]][["Times 30"]][[x]] <- linDecon_out(sim_sampled_dir_all$zero_5$`Times 30`[[x]],
                                                     sim_percent5)
}

for(x in paste("Rep", 1:5)){
  lin[["zero_50"]][["Times 30"]][[x]] <- linDecon_out(sim_sampled_dir_all$zero_50$`Times 30`[[x]],
                                                     sim_percent50)
}

#saveRDS(lin, file = "lin_zerofraction.rds")



