####Background####

# run time with sample size

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Large_samplesize")

load("/media/gcpeac/Myles/Deconvolution/AMP_scRNAseq_data_and_annotations.RData")

mk <- readRDS("../cellmarkers_new.rds") 
celseq_counts <- as.matrix(celseq_counts)

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

sim_out <- function(n){
  set.seed(3)
  sim_counts <- generate_samples(mk, n)
  sim_percent <- sim_counts / rowSums(sim_counts) * 100
  sim_pseudo <- simulate_bulk(mk, sim_counts)
  list(sim_counts = sim_counts,
       sim_percent = sim_percent,
       sim_pseudo)
}

large_sim <- list()

for(i in c(100, 500, 1000)){
  large_sim[[paste0(i, "_Samples")]] <- sim_out(i)
}

large_sim_sampled <- list()

for(i in c(100, 500, 1000)){
  for(x in paste("Rep", 1:5)){
    sim_counts <- large_sim[[paste0(i, "_Samples")]]$sim_counts
    large_sim_sampled[[paste("Times 30 and ", i, " Samples")]][[x]] <- 
      simulate_bulk(celseq_counts, 
                    sim_counts,
                    metadata$subclass,
                    times = 30)
  }
}

# add 5000 samples


large_sim[[paste0(5000, "_Samples")]] <- sim_out(5000)


# need to name pseudo

for(i in names(large_sim)){
  names(large_sim[[i]])[3] <- "sim_pseudo"
}

# get sim_sampled

for(x in paste("Rep", 1:5)){
  sim_counts <- large_sim[[paste0(5000, "_Samples")]]$sim_counts
  large_sim_sampled[[paste("Times 30 and ", 5000, " Samples")]][[x]] <- 
    simulate_bulk(celseq_counts, 
                  sim_counts,
                  metadata$subclass,
                  times = 30)
}

# save(large_sim, large_sim_sampled, file = "large_simulated_data.rdata")

####Cell Geometry####

library(cellGeometry)

geo_output <- function(n, rep){
  
  samples_sim <- large_sim_sampled[[paste("Times 30 and ", n, " Samples")]][[paste("Rep", rep)]]
  sim_counts <- large_sim[[paste0(n, "_Samples")]]$sim_counts
  sim_percent <- large_sim[[paste0(n, "_Samples")]]$sim_percent
  tic()
  out <- deconvolute(mk, samples_sim,
                     comp_amount = 1,
                     arith_mean = T,
                     count_space = T,
                     convert_bulk = F, use_filter = F, adjust_comp = T)
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

cellgeo <- list()

library(tictoc)

for(i in c(100, 500, 1000)){
  cellgeo[[paste("Times 30 and ", i, " Samples")]] <- list()
  for(x in paste("Rep", 1:5)){
    cellgeo[[paste("Times 30 and ", i, " Samples")]][[x]] <- geo_output(n = i, 
                                                                        rep = as.numeric(gsub("Rep ", "", x)))
  }
}

# addition of 5000 samples

for(x in paste("Rep", 1:5)){
  cellgeo[[paste("Times 30 and ", 5000, " Samples")]][[x]] <- geo_output(n = 5000, 
                                                                         rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(cellgeo, "cellgeo_output_large_new.rds")

####MuSiC####

library(Biobase)
library(MuSiC, lib.loc = "/usr/local/lib/R/site-library")
#MuSiC1


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


music_out <- function(n, rep){
  samples_sim <- large_sim_sampled[[paste("Times 30 and ", n, " Samples")]][[paste("Rep", rep)]]
  sim_percent <- large_sim[[paste0(n, "_Samples")]]$sim_percent
  
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

for(i in c(100, 500, 1000)){
  music1[[paste("Times 30 and ", i, " Samples")]] <- list()
  for(x in paste("Rep", 1:5)){
    music1[[paste("Times 30 and ", i, " Samples")]][[x]] <- music_out(n = i, 
                                                                      rep = as.numeric(gsub("Rep ", "", x)))
  }
}


# saveRDS(music1, "music1_output_large.rds")

# addition of 5000 samples

music1 <- readRDS("music1_output_large.rds")

for(x in paste("Rep", 1:4)){
  music1[[paste("Times 30 and ", 5000, " Samples")]][[x]] <- music_out(n = 5000, 
                                                                       rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(music1, "music1_output_large.rds")

# run 5th repeat

music1[[paste("Times 30 and ", 5000, " Samples")]][["Rep 5"]] <- music_out(n = 5000, 
                                                                           rep = 5)

# saveRDS(music1, "music1_output_large.rds")

####MuSiC2####

rm(list = setdiff(ls(), c("celseq_counts",
                          "large_sim",
                          "large_sim_sampled",
                          "metadata",
                          "mk")))

detach("package:MuSiC", unload = TRUE)
#restart
library(SingleCellExperiment)
library(MuSiC, lib.loc =  "/home/rachel/R/x86_64-pc-linux-gnu-library/4.2")
library(cellGeometry)

sce <- SingleCellExperiment(list(counts = celseq_counts[ ,is.na(metadata$subclass) == FALSE]),
                            colData = metadata[is.na(metadata$subclass) == FALSE, ])

music2_out <- function(n, rep){
  samples_sim <- large_sim_sampled[[paste("Times 30 and ", n, " Samples")]][[paste("Rep", rep)]]
  
  tic()
  out <- music_prop(bulk.mtx = samples_sim,
                    sc.sce = sce,
                    clusters = "subclass",
                    samples = "sample")
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  list(time = mk_time)
}

music2 <- list()
library(tictoc)
for(i in c(100, 500, 1000)){
  music2[[paste("Times 30 and ", i, " Samples")]] <- list()
  for(x in paste("Rep", 1:5)){
    music2[[paste("Times 30 and ", i, " Samples")]][[x]] <- music2_out(n = i, 
                                                                       rep = as.numeric(gsub("Rep ", "", x)))
  }
}


#saveRDS(music2, "music2_output_large.rds")

#addition of 5000 samples

music2 <- readRDS("music2_output_large.rds")

for(x in paste("Rep", 1:5)){
  music2[[paste("Times 30 and ", 5000, " Samples")]][[x]] <- music2_out(n = 5000, 
                                                                        rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(music2, "music2_output_large.rds")

####LinDeconSeq####

#restart R/workstation

library(LinDeconSeq)
library(cellGeometry)

load("large_simulated_data.rdata")

load("/media/gcpeac/Myles/Deconvolution/AMP_scRNAseq_data_and_annotations.RData")

celseq_counts <- as.matrix(celseq_counts)

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

dataSC <- celseq_counts[ , is.na(metadata$subclass) == FALSE]
meta <- metadata[is.na(metadata$subclass) == FALSE, ]

library(tidyr)
phes <- as.data.frame(meta %>% pivot_wider(subclass,
                                           names_from = cell_name,
                                           values_from = cell_name))
rownames(phes) <- phes$subclass
phes <- phes[ , -1]
phes <- as.matrix(phes)
phes <- ifelse(is.na(phes), 0, 1)

markerRes <- readRDS("../LinDeconSeq_signature.rds")

library(tictoc)
linDecon_out <- function(n, rep){
  samples_sim <- large_sim_sampled[[paste("Times 30 and ", n, " Samples")]][[paste("Rep", rep)]]
  sim_percent <- large_sim[[paste0(n, "_Samples")]]$sim_percent
  
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

linDecon <- list()

for(i in c(100, 500, 1000)){
  linDecon[[paste("Times 30 and ", i, " Samples")]] <- list()
  for(x in paste("Rep", 1:5)){
    linDecon[[paste("Times 30 and ", i, " Samples")]][[x]] <- linDecon_out(n = i, 
                                                                           rep = as.numeric(gsub("Rep ", "", x)))
  }
}


#saveRDS(linDecon, "LinDeconSeq_output_large.rds")

#addition of 5000 samples

linDecon <- readRDS("LinDeconSeq_output_large.rds")

for(x in paste("Rep", 1:5)){
  linDecon[[paste("Times 30 and ", 5000, " Samples")]][[x]] <- linDecon_out(n = 5000, 
                                                                            rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(linDecon, "LinDeconSeq_output_large.rds")

####DWLS####

library(DWLS)

load("Signature.rdata")
load("large_simulated_data.rdata")


DWLS_out <- function(n, rep){
  samples_sim <- large_sim_sampled[[paste("Times 30 and ", n, " Samples")]][[paste("Rep", rep)]]
  sim_percent <- large_sim[[paste0(n, "_Samples")]]$sim_percent
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
  
  # colnames(allCounts_DWLS) <- colnames(samples_sim)
  # colnames(allCounts_DWLS) <- colnames(samples_sim)
  # rownames(allCounts_DWLS) <- paste0("SC-", rownames(allCounts_DWLS))
  # propdf <- t(allCounts_DWLS) * 100
  # propdf <- propdf[ , match(colnames(sim_percent),
  #                           colnames(propdf))] #new addition
  # metric_percent <- metric_set(sim_percent, propdf)
  
  #hashed out for N = 5000
  
  list(
    # output = allCounts_DWLS,
    # metrics_percent = metric_percent,
    time = time)
  
}

library(tictoc)
library(cellGeometry)

DWLS <- list()

#100 samples

for(x in paste("Rep", 1:5)){
  DWLS[[paste("Times 30 and ", 100, " Samples")]][[x]] <- DWLS_out(n = 100,
                                                                   rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(DWLS, "DWLS_output_large.rds")



#500 samples

DWLS <- readRDS("DWLS_output_large.rds")
for(x in paste("Rep", 1:5)){
 DWLS[[paste("Times 30 and ", 500, " Samples")]][[x]] <- DWLS_out(n = 500,
                                                                  rep = as.numeric(gsub("Rep ", "", x)))
}

#saveRDS(DWLS, "DWLS_output_large.rds")

# 1000 samples

DWLS <- readRDS("DWLS_output_large.rds")

for(x in paste("Rep", 1:5)){
  DWLS[[paste("Times 30 and ", 1000, " Samples")]][[x]] <- DWLS_out(n = 1000,
                                                                   rep = as.numeric(gsub("Rep ", "", x)))
}

# saveRDS(DWLS, "DWLS_output_large.rds")

# 5000 samples

DWLS <- readRDS("DWLS_output_large.rds")

DWLS[[paste("Times 30 and ", 5000, " Samples")]][["Rep 1"]] <- DWLS_out(n = 5000,
                                                                    rep = 1)

#53115.568 sec elapsed
#saveRDS(DWLS, "DWLS_output_large.rds")

DWLS[[paste("Times 30 and ", 5000, " Samples")]][["Rep 2"]] <- DWLS_out(n = 5000,
                                                                        rep = 2)
#57636.192 sec elapsed

# saveRDS(DWLS, "DWLS_output_large.rds")

DWLS[[paste("Times 30 and ", 5000, " Samples")]][["Rep 3"]] <- DWLS_out(n = 5000,
                                                                        rep = 3)
#57265.903 sec elapsed
#saveRDS(DWLS, "DWLS_output_large.rds")

DWLS[[paste("Times 30 and ", 5000, " Samples")]][["Rep 4"]] <- DWLS_out(n = 5000,
                                                                        rep = 4)
#58439.155 sec elapsed
# saveRDS(DWLS, "DWLS_output_large.rds")

DWLS[[paste("Times 30 and ", 5000, " Samples")]][["Rep 5"]] <- DWLS_out(n = 5000,
                                                                        rep = 5)

saveRDS(DWLS, "DWLS_output_large.rds")






