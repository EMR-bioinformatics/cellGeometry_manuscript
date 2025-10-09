####Background####

# Using cellGeometry investigate different levels of noise
# This is different from the original script

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/Noise_updated")

sim_sampled_dir_all <- readRDS("../simulated_dirichlet_all.rds")
load("../simulated_dirichlet.rdata")

####Default noise generation####

noise_all <- list()
noise_all[["add"]] <- list()

for(i in paste("Rep", 1:5)){
  noise_all[["add"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]])
}

noise_all[["log"]] <- list()
for(i in paste("Rep", 1:5)){
  noise_all[["log"]][[i]] <- log_noise(sim_sampled_dir_all$`Times 30`[[i]])
}

noise_all[["sqrt"]] <- list()
for(i in paste("Rep", 1:5)){
  noise_all[["sqrt"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]])
}

noise_all[["shift"]] <- list()
for(i in paste("Rep", 1:5)){
  noise_all[["shift"]][[i]] <- shift_noise(sim_sampled_dir_all$`Times 30`[[i]])
}

#saveRDS(noise_all, file = "default_noise.rds")

noise_all <- readRDS("default_noise.rds")

####Default noise deconvolution####

mk <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/cellmarkers_new.rds")
mk_update <- updateMarkers(mk, nsubclass = 500, expfilter = 0.25)

library(cellGeometry)

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

cellgeo <- list()


library(tictoc)

for(i in names(noise_all)){
  cellgeo[[i]] <- list()

  for(x in names(noise_all[[i]])){
    cellgeo[[i]][[x]]<- geo_output(noise_all[[i]][[x]],
                                   mk_update)
  }
}

#saveRDS(cellgeo, file = "cellgeo_default_noise.rds")

####____different shift noise____####

shift_all <- list()
shift_all[["SD = 0.75"]] <- list()
for(i in paste("Rep", 1:5)){
  shift_all[["SD = 0.75"]][[i]] <- shift_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 0.75)
}

shift_all[["SD = 1"]] <- list()
for(i in paste("Rep", 1:5)){
  shift_all[["SD = 1"]][[i]] <- shift_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                               sd = 1)
}

shift_all[["SD = 1.5"]] <- list()
for(i in paste("Rep", 1:5)){
  shift_all[["SD = 1.5"]][[i]] <- shift_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                            sd = 1.5)
}

shift_all[["SD = 2"]] <- list()
for(i in paste("Rep", 1:5)){
  shift_all[["SD = 2"]][[i]] <- shift_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                              sd = 2)
}

#saveRDS(shift_all, "shift_noise.rds")

#####shift noise deconvolution#####

cellgeo_shift <- list()

for(i in names(shift_all)){
  cellgeo_shift[[i]] <- list()
  
  for(x in names(shift_all[[i]])){
    cellgeo_shift[[i]][[x]]<- geo_output(shift_all[[i]][[x]],
                                         mk_update)
  }
}

#saveRDS(cellgeo_shift, "cellgeo_shift_noise.rds")

#####shift noise deconvolution - npass = 2#####

geo_outputv2 <- function(mat, mk){
  
  samples_sim <- mat
  tic()
  out <- deconvolute(mk, samples_sim,
                     arith_mean = T,
                     count_space = T, convert_bulk = F, 
                     use_filter = FALSE,
                     weight_method = "equal",
                     var_cutoff = 2.5, 
                     npass = 2)
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

cellgeo_shift_npass <- list()

cellgeo_shift_npass[["SD = 0"]] <- list()

for(i in names(sim_sampled_dir_all$`Times 30`)){
  cellgeo_shift_npass[["SD = 0"]][[i]] <- geo_outputv2(sim_sampled_dir_all$`Times 30`[[i]],
                                                       mk_update)
}

for(i in names(noise_all[["shift"]])){
  cellgeo_shift_npass[["SD = 0.5"]][[i]] <- geo_outputv2(noise_all[["shift"]][[i]],
                                                         mk_update)
}


for(i in names(shift_all)){
  cellgeo_shift_npass[[i]] <- list()
  
  for(x in names(shift_all[[i]])){
    cellgeo_shift_npass[[i]][[x]]<- geo_outputv2(shift_all[[i]][[x]],
                                         mk_update)
  }
}

#saveRDS(cellgeo_shift_npass, "cellgeo_shift_npass_var2point5.rds")


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

DWLS_noise <- list()

DWLS_noise[["SD = 0"]] <- list()

for(i in names(noise_all[["shift"]])){
  DWLS_noise[["SD = 0.5"]][[i]] <- DWLS_out(noise_all[["shift"]][[i]])
}

for(x in names(shift_all)){
  DWLS_noise[[x]] <- list()
  
  for(i in names(shift_all[[x]])){
    DWLS_noise[[x]][[i]] <- DWLS_out(shift_all[[x]][[i]])
  }
}

#saveRDS(DWLS_noise, "DWLS_shift_noise.rds")


####____different add noise____####

add_all <- list()
add_all[["SD = 250"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 250"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                               sd = 250)
}

add_all[["SD = 500"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 500"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                            sd = 500)
}

add_all[["SD = 750"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 750"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                              sd = 750)
}

add_all[["SD = 1000"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 1000"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                            sd = 1000)
}

#saveRDS(add_all, "add_noise.rds")

#1000 too much- 500 maximum 

add_all <- readRDS("add_noise.rds")

add_all[["SD = 200"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 200"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 200)
}

add_all[["SD = 300"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 300"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 300)
}

add_all[["SD = 400"]] <- list()
for(i in paste("Rep", 1:5)){
  add_all[["SD = 400"]][[i]] <- add_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 400)
}

#saveRDS(add_all, "add_noise.rds")

#####add noise deconvolution####

cellgeo_add <- list()

for(i in names(add_all)){
  cellgeo_add[[i]] <- list()
  
  for(x in names(add_all[[i]])){
    cellgeo_add[[i]][[x]]<- geo_output(add_all[[i]][[x]],
                                         mk_update)
  }
}

#saveRDS(cellgeo_add, "cellgeo_add_noise.rds")

for(i in names(add_all)[5:7]){
  cellgeo_add[[i]] <- list()
  
  for(x in names(add_all[[i]])){
    cellgeo_add[[i]][[x]]<- geo_output(add_all[[i]][[x]],
                                       mk_update)
  }
}

#saveRDS(cellgeo_add, "cellgeo_add_noise.rds")

#####DWLS####

DWLS_add <- list()

DWLS_add[["SD = 0"]] <- list()

for(i in names(noise_all[["add"]])){
  DWLS_add[["SD = 100"]][[i]] <- DWLS_out(noise_all[["add"]][[i]])
}

for(x in names(add_all)){
  DWLS_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    DWLS_add[[x]][[i]] <- DWLS_out(add_all[[x]][[i]])
  }
}

#saveRDS(DWLS_add, "DWLS_add_noise.rds")

for(x in names(add_all)[5:7]){
  DWLS_add[[x]] <- list()
  
  for(i in names(add_all[[x]])){
    DWLS_add[[x]][[i]] <- DWLS_out(add_all[[x]][[i]])
  }
}

#saveRDS(DWLS_add, "DWLS_add_noise.rds")

####____different log noise____####

log_all <- list()
log_all[["SD = 0.25"]] <- list()
for(i in paste("Rep", 1:5)){
  log_all[["SD = 0.25"]][[i]] <- log_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 0.25)
}

log_all[["SD = 0.5"]] <- list()
for(i in paste("Rep", 1:5)){
  log_all[["SD = 0.5"]][[i]] <- log_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 0.5)
}

log_all[["SD = 0.75"]] <- list()
for(i in paste("Rep", 1:5)){
  log_all[["SD = 0.75"]][[i]] <- log_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 0.75)
}

log_all[["SD = 1"]] <- list()
for(i in paste("Rep", 1:5)){
  log_all[["SD = 1"]][[i]] <- log_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 1)
}

#saveRDS(log_all, "log_noise.rds")

#####log noise deconvolution####

cellgeo_log <- list()

for(i in names(log_all)){
  cellgeo_log[[i]] <- list()
  
  for(x in names(log_all[[i]])){
    cellgeo_log[[i]][[x]]<- geo_output(log_all[[i]][[x]],
                                       mk_update)
  }
}

#saveRDS(cellgeo_log, "cellgeo_log_noise.rds")

#####DWLS####

DWLS_log <- list()

DWLS_log[["SD = 0"]] <- list()

for(i in names(noise_all[["log"]])){
  DWLS_log[["SD = 0.1"]][[i]] <- DWLS_out(noise_all[["log"]][[i]])
}

for(x in names(log_all)){
  DWLS_log[[x]] <- list()
  
  for(i in names(log_all[[x]])){
    DWLS_log[[x]][[i]] <- DWLS_out(log_all[[x]][[i]])
  }
}

#saveRDS(DWLS_log, "DWLS_log_noise.rds")

####____different sqrt noise____####
sqrt_all <- list()
sqrt_all[["SD = 250"]] <- list()
for(i in paste("Rep", 1:5)){
  sqrt_all[["SD = 250"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 250)
}

sqrt_all[["SD = 500"]] <- list()
for(i in paste("Rep", 1:5)){
  sqrt_all[["SD = 500"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 500)
}

sqrt_all[["SD = 750"]] <- list()
for(i in paste("Rep", 1:5)){
  sqrt_all[["SD = 750"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                          sd = 750)
}

sqrt_all[["SD = 1000"]] <- list()
for(i in paste("Rep", 1:5)){
  sqrt_all[["SD = 1000"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 1000)
}

#saveRDS(sqrt_all, "sqrt_noise.rds")

#1000 too much. 100 maximum 

sqrt_all[["SD = 10"]] <- list()
for(i in paste("Rep", 1:5)){
  sqrt_all[["SD = 10"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                             sd = 10)
}

sqrt_all[["SD = 25"]] <- list()
for(i in paste("Rep", 1:5)){
  sqrt_all[["SD = 25"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 25)
}

sqrt_all[["SD = 50"]] <- list()
for(i in paste("Rep", 1:5)){
  sqrt_all[["SD = 50"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 50)
}

sqrt_all[["SD = 75"]] <- list()
for(i in paste("Rep", 1:5)){
  sqrt_all[["SD = 75"]][[i]] <- sqrt_noise(sim_sampled_dir_all$`Times 30`[[i]],
                                           sd = 75)
}

#saveRDS(sqrt_all, "sqrt_noise.rds")

#####sqrt noise deconvolution####

cellgeo_sqrt <- list()

for(i in names(sqrt_all)){
  cellgeo_sqrt[[i]] <- list()
  
  for(x in names(sqrt_all[[i]])){
    cellgeo_sqrt[[i]][[x]]<- geo_output(sqrt_all[[i]][[x]],
                                       mk_update)
  }
}

#saveRDS(cellgeo_sqrt, "cellgeo_sqrt_noise.rds")

for(i in names(sqrt_all)[5:8]){
  cellgeo_sqrt[[i]] <- list()
  
  for(x in names(sqrt_all[[i]])){
    cellgeo_sqrt[[i]][[x]]<- geo_output(sqrt_all[[i]][[x]],
                                        mk_update)
  }
}

#saveRDS(cellgeo_sqrt, "cellgeo_sqrt_noise.rds")

#####DWLS####

DWLS_sqrt <- list()

DWLS_sqrt[["SD = 0"]] <- list()

for(i in names(noise_all[["sqrt"]])){
  DWLS_sqrt[["SD = 100"]][[i]] <- DWLS_out(noise_all[["sqrt"]][[i]])
}

for(x in names(sqrt_all)){
  DWLS_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    DWLS_sqrt[[x]][[i]] <- DWLS_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(DWLS_sqrt, "DWLS_sqrt_noise.rds")

for(x in names(sqrt_all)[5:8]){
  DWLS_sqrt[[x]] <- list()
  
  for(i in names(sqrt_all[[x]])){
    DWLS_sqrt[[x]][[i]] <- DWLS_out(sqrt_all[[x]][[i]])
  }
}

#saveRDS(DWLS_sqrt, "DWLS_sqrt_noise.rds")

