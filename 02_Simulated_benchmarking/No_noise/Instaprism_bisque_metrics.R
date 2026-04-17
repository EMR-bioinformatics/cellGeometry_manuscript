#retrieving RMSE, Rsq and CCC for instaprism and bisque 

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/metrics")

library(cellGeometry)
library(InstaPrism)
library(BisqueRNA)

rmse <- function(obs, pred) {
  sqrt(mean((pred - obs)^2))
}

Rsq <- function(obs, pred) {
  rss <- sum((pred - obs)^2)
  tss <- sum((obs - mean(obs))^2)
  1 - rss/tss
}

metric_all <- function(df, sim, method, rep){
  sub <- df[df$Method == method & df$Rep == rep, ]
  stopifnot(all(paste(sub$Patient, sub$Cluster) == 
                  paste(sim$Patient, sim$Cluster)))

  data.frame("RMSE" = rmse(sim$Percentage, sub$Percentage),
             "Rsq" = Rsq(sim$Percentage, sub$Percentage),
             "Method" = method,
             "Rep" = rep)
}

library(DescTools)

CCC_all <- function(df, sim, method, rep){
  sub <- df[df$Method == method & df$Rep == rep, ]
  
  stopifnot(all(paste(sub$Patient, sub$Cluster) == 
                  paste(sim$Patient, sim$Cluster)))
  
  CCC(sub$Percentage, sim$Percentage)
}

####tabula####

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/tabula_simulated_workstation2.rdata")

sim_percent_long <- as.data.frame(sim_percent)

sim_percent_long$Patient <- rownames(sim_percent_long)

sim_percent_long <- as.data.frame(pivot_longer(sim_percent_long,
                                               cols = -Patient,
                                               names_to = "Cluster",
                                               values_to = "Percentage"))

insta <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/instaprism_tabula/instaprism_tabula_output.rds")
bisque <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/Bisque_tabula/bisque_tabula_output.rds")

insta_long <- lapply(names(insta$`Times 3`), function(x){
  df <- as.data.frame(t(insta$`Times 3`[[x]]$output@Post.ini.ct@theta) * 100)
  df <- df[ , match(colnames(sim_percent),
                    colnames(df))]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "InstaPrism"
  test$Rep <- x
  test
})

celltypedf <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/DWLS_celltype_convert.rds")
##subset dataset same as DWLS

bisque_long <- lapply(names(bisque$`Times 3`), function(x){
  df <- as.data.frame(t(bisque[["Times 3"]][[x]][["output"]][["bulk.props"]]) * 100)
  colnames(df) <- celltypedf$Orig[match(colnames(df),
                                        celltypedf$convert)]
  df <- df[ , match(colnames(sim_percent),
                    colnames(df))]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "Bisque"
  test$Rep <- x
  test
  
})

comb <- rbind(do.call(rbind, insta_long),
              do.call(rbind, bisque_long))

CCC_comb <- list()

for(i in c("InstaPrism", "Bisque")){
  for(x in paste("Rep", 1:5)){
    CCC_comb[[paste0(i, "_", gsub(" ", "", x))]] <- CCC_all(comb, sim_percent_long,
                                                            i, x)
  }
}

#saveRDS(CCC_comb, file = "insta_bisque_tabula_ccc.rds")

CCC_comb_res <- data.frame()

for(i in names(CCC_comb)){
  temp <- CCC_comb[[i]]$rho.c
  temp$Method <- gsub("_.*", "", i)
  temp$Rep <- gsub(".*_", "", i)
  temp$Rep <- gsub("Rep", "", temp$Rep)
  CCC_comb_res <- rbind(CCC_comb_res,
                        temp)
}

tabula_comb <- list()

for(i in c("InstaPrism", "Bisque")){
  for(x in paste("Rep", 1:5)){
    tabula_comb[[paste0(i, "_", gsub(" ", "", x))]] <- metric_all(comb, sim_percent_long,
                                                                  i, x)
  }
}

tabula_comb <- do.call(rbind, tabula_comb)
tabula_comb$Data <- "Tabula Sapiens"

tabula_comb$CCC <- CCC_comb_res$est[match(paste0(tabula_comb$Method, gsub("Rep ", "", tabula_comb$Rep)),
                                          paste0(CCC_comb_res$Method, CCC_comb_res$Rep))]

####Brain####

rm(list = setdiff(ls(), c("CCC_all", "metric_all",
                          "rmse", "Rsq",
                          "tabula_comb")))

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/brain_simulated_merged.rdata")

insta <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/Instaprism_brain/instaprism_neuron_output.rds")
insta_NN <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/Instaprism_brain/instaprism_NN_output.rds")

bisque <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/Bisque_brain/bisque_neuron_output.rds")
bisque_NN <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/Bisque_brain/bisque_NN_output.rds")

sim_percent_neuron <- sim_percent_merge[ , grepl("Human", colnames(sim_percent_merge))]
sim_percent_NN <- sim_percent_merge[ , !grepl("Human", colnames(sim_percent_merge))]

insta_long <- lapply(names(insta$`Times 1`), function(x){
  df <- as.data.frame(rbind(t(insta$`Times 1`[[x]]$output@Post.ini.ct@theta) * 100))
  df <- df[ , match(colnames(sim_percent_neuron),
                    colnames(df))]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "InstaPrism"
  test$Rep <- x
  test
})

insta_long_NN <- lapply(names(insta_NN$`Times 1`), function(x){
  df <- as.data.frame(rbind(t(insta_NN$`Times 1`[[x]]$output@Post.ini.ct@theta) * 100))
  df <- df[ , match(colnames(sim_percent_NN),
                    colnames(df))]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "InstaPrism"
  test$Rep <- x
  test
})

celltype <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/DWLS_neuron_convert.rds")
#subset dataset same as DWLS

bisque_long <- lapply(names(bisque$`Times 1`), function(x){
  df <- as.data.frame(t(bisque[["Times 1"]][[x]][["output"]][["bulk.props"]]) * 100)
  colnames(df) <- celltype$Orig[match(colnames(df),
                                        celltype$convert)]
  df <- df[ , match(colnames(sim_percent_neuron),
                    colnames(df))]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "Bisque"
  test$Rep <- x
  test
  
})

bisque_long_NN <- lapply(names(bisque_NN$`Times 1`), function(x){
  df <- as.data.frame(t(bisque_NN[["Times 1"]][[x]][["output"]][["bulk.props"]]) * 100)
  df <- df[ , match(colnames(sim_percent_NN),
                    colnames(df))]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "Bisque"
  test$Rep <- x
  test
  
})

comb <- rbind(do.call(rbind, insta_long),
              do.call(rbind, insta_long_NN),
              do.call(rbind, bisque_long),
              do.call(rbind, bisque_long_NN))

comb$Data <- ifelse(grepl("Human", comb$Cluster), 
                    "Human Brain Cell Atlas (Neuron)", 
                    "Human Brain Cell Atlas (Non-neuronal)")

sim_percent_long <- as.data.frame(sim_percent_merge)

sim_percent_long$Patient <- rownames(sim_percent_long)

sim_percent_long <- as.data.frame(pivot_longer(sim_percent_long,
                                               cols = -Patient,
                                               names_to = "Cluster",
                                               values_to = "Percentage"))

sim_percent_long_neuron <- sim_percent_long[grepl("Human", sim_percent_long$Cluster), ]
sim_percent_long_NN <- sim_percent_long[!grepl("Human", sim_percent_long$Cluster), ]

CCC_comb_neuron <- list()

for(i in c("InstaPrism", "Bisque")){
  for(x in paste("Rep", 1:3)){
    CCC_comb_neuron[[paste0(i, "_", gsub(" ", "", x))]] <- CCC_all(comb[grepl("Human", comb$Cluster), ],
                                                                       sim_percent_long_neuron,
                                                                       i, x)
  }
}

#saveRDS(CCC_comb_neuron, file = "insta_bisque_neuron_ccc.rds")

neuron_comb <- list()

for(i in c("InstaPrism", "Bisque")){
  for(x in paste("Rep", 1:3)){
    neuron_comb[[paste0(i, "_", gsub(" ", "", x))]] <- metric_all(comb[grepl("Human", comb$Cluster), ],
                                                                  sim_percent_long_neuron,
                                                                  i, x)
  }
}

neuron_comb <- do.call(rbind, neuron_comb)
neuron_comb$Data <- "Human Brain Cell Atlas (Neuron)"

CCC_comb_res_neuron <- data.frame()

for(i in names(CCC_comb_neuron)){
  temp <- CCC_comb_neuron[[i]]$rho.c
  temp$Method <- gsub("_.*", "", i)
  temp$Rep <- gsub(".*_", "", i)
  temp$Rep <- gsub("Rep", "", temp$Rep)
  CCC_comb_res_neuron <- rbind(CCC_comb_res_neuron,
                               temp)
}

neuron_comb$CCC <- CCC_comb_res_neuron$est[match(paste0(neuron_comb$Method, gsub("Rep ", "", neuron_comb$Rep)),
                                                 paste0(CCC_comb_res_neuron$Method, CCC_comb_res_neuron$Rep))]

CCC_comb_NN <- list()

for(i in c("InstaPrism", "Bisque")){
  for(x in paste("Rep", 1:3)){
    CCC_comb_NN[[paste0(i, "_", gsub(" ", "", x))]] <- CCC_all(comb[!grepl("Human", comb$Cluster), ],
                                                                   sim_percent_long_NN,
                                                                   i, x)
  }
}

#saveRDS(CCC_comb_NN, file = "insta_bisque_NN_ccc.rds")

NN_comb <- list()

for(i in c("InstaPrism", "Bisque")){
  for(x in paste("Rep", 1:3)){
    NN_comb[[paste0(i, "_", gsub(" ", "", x))]] <- metric_all(comb[!grepl("Human", comb$Cluster), ],
                                                              sim_percent_long_NN,
                                                              i, x)
  }
}

NN_comb <- do.call(rbind, NN_comb)
NN_comb$Data <- "Human Brain Cell Atlas (Non-neuronal)"

CCC_comb_res_NN <- data.frame()

for(i in names(CCC_comb_NN)){
  temp <- CCC_comb_NN[[i]]$rho.c
  temp$Method <- gsub("_.*", "", i)
  temp$Rep <- gsub(".*_", "", i)
  temp$Rep <- gsub("Rep", "", temp$Rep)
  CCC_comb_res_NN <- rbind(CCC_comb_res_NN,
                           temp)
}

NN_comb$CCC <- CCC_comb_res_NN$est[match(paste0(NN_comb$Method, gsub("Rep ", "", NN_comb$Rep)),
                                                 paste0(CCC_comb_res_NN$Method, CCC_comb_res_NN$Rep))]

####combine####

rm(list = setdiff(ls(), c("CCC_all", "metric_all",
                          "rmse", "Rsq",
                          "tabula_comb", "neuron_comb",
                          "NN_comb")))


comb <- rbind(tabula_comb,
              neuron_comb,
              NN_comb)

comb$Data <- factor(comb$Data,
                    levels = c("Tabula Sapiens",
                               "Human Brain Cell Atlas (Neuron)",
                               "Human Brain Cell Atlas (Non-neuronal)"))
                      
comb_mean <- comb %>%
  group_by(Data, Method) %>%
  summarise("RMSE_mean" = mean(RMSE),
            "Rsq_mean" = mean(Rsq),
            "CCC_mean" = mean(CCC),
            "N" = n(),
            "RMSE_SD" = sd(RMSE),
            "Rsq_SD" = sd(Rsq),
            "CCC_SD" = sd(CCC),
            "RMSE_SEM" = RMSE_SD/sqrt(N),
            "Rsq_SEM" = Rsq_SD/sqrt(N),
            "CCC_SEM" = CCC_SD/sqrt(N))

comb_mean$RMSE <- paste(formatC(comb_mean$RMSE_mean, 3),
                        "±", 
                        formatC(comb_mean$RMSE_SEM, 3))

comb_mean$Rsq <- paste(formatC(comb_mean$Rsq_mean, 3),
                       "±", 
                       formatC(comb_mean$Rsq_SEM, 3))

comb_mean$CCC <- paste(formatC(comb_mean$CCC_mean, 3),
                        "±", 
                        formatC(comb_mean$CCC_SEM, 3))

# write.table(comb_mean[ , c("Data", "Method", "RMSE", "Rsq", "CCC")],
#             file = "insta_bisque_metrics.txt",
#             append = FALSE, sep = "\t",
#             col.names = NA, row.names = TRUE)

