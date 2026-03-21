#calculate CCC and visualise Bland-Altman plots

####AMP####

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet")

library(cellGeometry)

cellgeo <- readRDS("cellgeo_dirichlet_output500.rds")

music2 <- readRDS("Music_dirichlet_output.rds")
DWLS <- readRDS("DWLS_dirichlet_output.rds")
Lin <- readRDS("Lin_dirichlet_output.rds")

load("simulated_dirichlet.rdata")

arrange_data <- function(df, method, rep){
  
  colnames(df) <- gsub("\\.", "-", colnames(df))
  df <- df[ , colnames(sim_percent)]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- method
  test$Rep <- rep
  test
}

library(tidyverse)

cellgeo_percent <- lapply(names(cellgeo$`Times 30`), function(x){
  arrange_data(as.data.frame(cellgeo$`Times 30`[[x]]$output$subclass$percent), "CellGeometry", x)
})

music_percent <- lapply(names(music2$`Times 30`), function(x){
  arrange_data(as.data.frame(music2$`Times 30`[[x]]$output$Est.prop.weighted * 100), "MuSiC", x)
})

DWLS_percent <- lapply(names(DWLS$`Times 30`), function(x){
  arrange_data(as.data.frame(t(DWLS$`Times 30`[[x]]$output  * 100)), "DWLS", x)
})

Lin_percent <- lapply(names(Lin$`Times 30`), function(x){
  arrange_data(as.data.frame(Lin$`Times 30`[[x]]$output * 100), "LinDeconSeq", x)
})

sim_percent_long <- as.data.frame(sim_percent)

sim_percent_long$Patient <- rownames(sim_percent_long)

sim_percent_long <- as.data.frame(pivot_longer(sim_percent_long,
                                               cols = -Patient,
                                               names_to = "Cluster",
                                               values_to = "Percentage"))

comb <- rbind(do.call(rbind, cellgeo_percent),
              do.call(rbind, music_percent),
              do.call(rbind, DWLS_percent),
              do.call(rbind, Lin_percent))

all(paste(comb$Patient[comb$Method == "CellGeometry" & comb$Rep == "Rep 1"], 
          comb$Cluster[comb$Method == "CellGeometry" & comb$Rep == "Rep 1"]) == 
      paste(sim_percent_long$Patient, sim_percent_long$Cluster))
#TRUE

library(DescTools)

CCC(comb$Percentage[comb$Method == "CellGeometry" & comb$Rep == "Rep 1"], sim_percent_long$Percentage)
# 0.973332 

CCC_all <- function(method, rep){
  sub <- comb[comb$Method == method & comb$Rep == rep, ]
  
  stopifnot(all(paste(sub$Patient, sub$Cluster) == 
                  paste(sim_percent_long$Patient, sim_percent_long$Cluster)))
  
  CCC(sub$Percentage, sim_percent_long$Percentage)
}

CCC_comb <- list()

for(i in c("CellGeometry", "MuSiC", "DWLS", "LinDeconSeq")){
  for(x in paste("Rep", 1:5)){
    CCC_comb[[paste0(i, "_", gsub(" ", "", x))]] <- CCC_all(i, x)
  }
}

#saveRDS(CCC_comb, file = "AMP_CCC_out.rds")

CCC_comb_res <- data.frame()

for(i in names(CCC_comb)){
  temp <- CCC_comb[[i]]$rho.c
  temp$Method <- gsub("_.*", "", i)
  temp$Rep <- gsub(".*_", "", i)
  temp$Rep <- gsub("Rep", "", temp$Rep)
  CCC_comb_res <- rbind(CCC_comb_res,
                        temp)
}

CCC_comb_res$Method <- factor(CCC_comb_res$Method,
                              levels = c("CellGeometry",
                                         "MuSiC",
                                         "DWLS",
                                         "LinDeconSeq"))

CCC_comb_mean <- CCC_comb_res %>%
  group_by(Method) %>%
  summarise("Mean" = mean(est),
            "SD" = sd(est),
            "N" = n(),
            "SEM" = SD/sqrt(N))

CCC_comb_mean$Mean <- formatC(CCC_comb_mean$Mean, 3)
CCC_comb_mean$SEM <- formatC(CCC_comb_mean$SEM, 3)

# write.table(CCC_comb_mean[ , c("Method", "Mean", "SEM")],
#             file = "CCC/CCC_comb_res_mean.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

####Typist####

rm(list = ls())

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Cell_typist/workstation2/Dirichlet")

cellgeo <- readRDS("cellgeo_dirichlet_output.rds")
music2 <- readRDS("music2_output.rds")
DWLS <- readRDS("DWLS_dirichlet_output.rds")
Lin <- readRDS("LinDeconSeq_output.rds")

load("simulated_dirichlet.rdata")

arrange_data <- function(df, method, rep){
  
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
  
  missing <- setdiff(colnames(sim_percent), colnames(df))
  
  if(length(missing) > 0){
    for(i in missing){
      df[ , i] <- 0
    }
  }
  
  
  df <- df[ , colnames(sim_percent)]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- method
  test$Rep <- rep
  test
}

cellgeo_percent <- lapply(names(cellgeo$`Times 30`), function(x){
  arrange_data(as.data.frame(cellgeo$`Times 30`[[x]]$output$subclass$percent), "CellGeometry", x)
})

music_percent <- lapply(names(music2$`Times 30`), function(x){
  arrange_data(as.data.frame(music2$`Times 30`[[x]]$output$Est.prop.weighted * 100), "MuSiC", x)
})

DWLS_percent <- lapply(names(DWLS$`Times 30`), function(x){
  arrange_data(as.data.frame(t(DWLS$`Times 30`[[x]]$output  * 100)), "DWLS", x)
})

Lin_percent <- lapply(names(Lin$`Times 30`), function(x){
  arrange_data(as.data.frame(Lin$`Times 30`[[x]]$output * 100), "LinDeconSeq", x)
})

sim_percent_long <- as.data.frame(sim_percent)

sim_percent_long$Patient <- rownames(sim_percent_long)

sim_percent_long <- as.data.frame(pivot_longer(sim_percent_long,
                                               cols = -Patient,
                                               names_to = "Cluster",
                                               values_to = "Percentage"))

comb <- rbind(do.call(rbind, cellgeo_percent),
              do.call(rbind, music_percent),
              do.call(rbind, DWLS_percent),
              do.call(rbind, Lin_percent))

CCC_all <- function(method, rep){
  sub <- comb[comb$Method == method & comb$Rep == rep, ]
  
  stopifnot(all(paste(sub$Patient, sub$Cluster) == 
                  paste(sim_percent_long$Patient, sim_percent_long$Cluster)))
  
  CCC(sub$Percentage, sim_percent_long$Percentage)
}

CCC_comb <- list()

for(i in c("CellGeometry", "MuSiC", "DWLS", "LinDeconSeq")){
  for(x in paste("Rep", 1:5)){
    CCC_comb[[paste0(i, "_", gsub(" ", "", x))]] <- CCC_all(i, x)
  }
}

#saveRDS(CCC_comb, file = "typist_CCC_out.rds")

CCC_comb_res <- data.frame()

for(i in names(CCC_comb)){
  temp <- CCC_comb[[i]]$rho.c
  temp$Method <- gsub("_.*", "", i)
  temp$Rep <- gsub(".*_", "", i)
  temp$Rep <- gsub("Rep", "", temp$Rep)
  CCC_comb_res <- rbind(CCC_comb_res,
                        temp)
}

CCC_comb_res$Method <- factor(CCC_comb_res$Method,
                              levels = c("CellGeometry",
                                         "MuSiC",
                                         "DWLS",
                                         "LinDeconSeq"))

CCC_comb_mean <- CCC_comb_res %>%
  group_by(Method) %>%
  summarise("Mean" = mean(est),
            "SD" = sd(est),
            "N" = n(),
            "SEM" = SD/sqrt(N))

CCC_comb_mean$Mean <- formatC(CCC_comb_mean$Mean, 3)
CCC_comb_mean$SEM <- formatC(CCC_comb_mean$SEM, 3)

# write.table(CCC_comb_mean[ , c("Method", "Mean", "SEM")],
#             file = "CCC/CCC_comb_res_mean.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

####Tabula####

rm(list = ls())

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2")

cellgeo <- readRDS("cellgeo_dirichlet_output_500.rds")
music2 <- readRDS("music_dirichlet_output.rds")
DWLS <- readRDS("DWLS_dirichlet_output.rds")

celltypedf <- readRDS("DWLS_celltype_convert.rds")

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/tabula_simulated_workstation2.rdata")

arrange_data <- function(df, method, rep){
  
  common <- intersect(colnames(sim_percent), colnames(df))
  df <- df[ , colnames(df) %in% common]
  
  missing <- setdiff(colnames(sim_percent), colnames(df))
  
  if(length(missing) > 0){
    for(i in missing){
      df[ , i] <- 0
    }
  }
  
  
  df <- df[ , colnames(sim_percent)]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- method
  test$Rep <- rep
  test
}

cellgeo_percent <- lapply(names(cellgeo$`Times 3`), function(x){
  arrange_data(as.data.frame(cellgeo$`Times 3`[[x]]$output$subclass$percent), "CellGeometry", x)
})

music_percent <- lapply(names(music2$`Times 3`), function(x){
  arrange_data(as.data.frame(music2$`Times 3`[[x]]$output$Est.prop.weighted * 100), "MuSiC", x)
})

arrange_DWLS <- function(df, rep){
  
  df <- as.data.frame(t(df) * 100)
  
  colnames(df) <- celltypedf$Orig[match(colnames(df),
                                        celltypedf$convert)]
  
  common <- intersect(colnames(sim_percent), colnames(df))
  df <- df[ , colnames(df) %in% common]
  
  missing <- setdiff(colnames(sim_percent), colnames(df))
  
  if(length(missing) > 0){
    for(i in missing){
      df[ , i] <- 0
    }
  }
  
  
  df <- df[ , colnames(sim_percent)]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "DWLS"
  test$Rep <- rep
  test
  
}


DWLS_percent <- lapply(names(DWLS$`Times 3`), function(x){
  arrange_DWLS(as.data.frame(DWLS$`Times 3`[[x]]$output), x)
})

sim_percent_long <- as.data.frame(sim_percent)

sim_percent_long$Patient <- rownames(sim_percent_long)

sim_percent_long <- as.data.frame(pivot_longer(sim_percent_long,
                                               cols = -Patient,
                                               names_to = "Cluster",
                                               values_to = "Percentage"))

comb <- rbind(do.call(rbind, cellgeo_percent),
              do.call(rbind, music_percent),
              do.call(rbind, DWLS_percent))

CCC_all <- function(method, rep){
  sub <- comb[comb$Method == method & comb$Rep == rep, ]
  
  stopifnot(all(paste(sub$Patient, sub$Cluster) == 
                  paste(sim_percent_long$Patient, sim_percent_long$Cluster)))
  
  CCC(sub$Percentage, sim_percent_long$Percentage)
}

CCC_comb <- list()

for(i in c("CellGeometry", "MuSiC", "DWLS")){
  for(x in paste("Rep", 1:5)){
    CCC_comb[[paste0(i, "_", gsub(" ", "", x))]] <- CCC_all(i, x)
  }
}

#saveRDS(CCC_comb, file = "Tabula_CCC_out.rds")

CCC_comb_res <- data.frame()

for(i in names(CCC_comb)){
  temp <- CCC_comb[[i]]$rho.c
  temp$Method <- gsub("_.*", "", i)
  temp$Rep <- gsub(".*_", "", i)
  temp$Rep <- gsub("Rep", "", temp$Rep)
  CCC_comb_res <- rbind(CCC_comb_res,
                        temp)
}

CCC_comb_res$Method <- factor(CCC_comb_res$Method,
                              levels = c("CellGeometry",
                                         "MuSiC",
                                         "DWLS"))

CCC_comb_mean <- CCC_comb_res %>%
  group_by(Method) %>%
  summarise("Mean" = mean(est),
            "SD" = sd(est),
            "N" = n(),
            "SEM" = SD/sqrt(N))

CCC_comb_mean$Mean <- formatC(CCC_comb_mean$Mean, 3)
CCC_comb_mean$SEM <- formatC(CCC_comb_mean$SEM, 3)

# write.table(CCC_comb_mean[ , c("Method", "Mean", "SEM")],
#             file = "CCC/CCC_comb_res_mean.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

####Brain####

rm(list = ls())

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain")

library(cellGeometry)

cellgeo <- readRDS("cellgeo_dirichlet_output.rds")
music2 <- readRDS("music_dirichlet_output.rds")
music2_NN <- readRDS("music_dirichlet_output_NN.rds")
DWLS_NN <- readRDS("DWLS_dirichlet_output_NN.rds")
Lin_NN <- readRDS("Lin_dirichlet_output_NN.rds")

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/brain_simulated_merged.rdata")

cellgeo_percent <- lapply(paste("Rep", 1:3), function(x){
  df <- as.data.frame(cellgeo$`Times 1`[[x]]$output$subclass$percent)
  
  df <- df[ , colnames(sim_percent_merge)]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "CellGeometry"
  test$Rep <- x
  test
  
})

music_percent <- lapply(paste("Rep", 1:3), function(x){
  df <- cbind(as.data.frame(music2$`Times 1`[[x]]$output$Est.prop.weighted * 100),
              as.data.frame(music2_NN$`Times 1`[[x]]$output$Est.prop.weighted * 100))
  
  common <- intersect(colnames(sim_percent_merge), colnames(df))
  df <- df[ , colnames(df) %in% common]
  
  missing <- setdiff(colnames(sim_percent_merge), colnames(df))
  
  if(length(missing) > 0){
    for(i in missing){
      df[ , i] <- 0
    }
  }
  
  
  df <- df[ , colnames(sim_percent_merge)]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- "MuSiC"
  test$Rep <- x
  test
})

sim_percent_neuron <- sim_percent_merge[ , grepl("Human", colnames(sim_percent_merge))]
sim_percent_NN <- sim_percent_merge[ , !grepl("Human", colnames(sim_percent_merge))]

arrange_NN <- function(df, method, rep){
  colnames(df) <- gsub("_", " ", colnames(df))
  df <- df[ , colnames(sim_percent_NN)]
  
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- method
  test$Rep <- rep
  test
}

DWLS_percent_NN <- lapply(names(DWLS_NN$`Times 1`), function(x){
  arrange_NN(as.data.frame(t(DWLS_NN$`Times 1`[[x]]$output  * 100)), "DWLS", x)
})

Lin_percent_NN <- lapply(names(Lin_NN$`Times 1`), function(x){
  arrange_NN(as.data.frame(Lin_NN$`Times 1`[[x]]$output * 100), "LinDeconSeq", x)
})

sim_percent_long <- as.data.frame(sim_percent_merge)

sim_percent_long$Patient <- rownames(sim_percent_long)

sim_percent_long <- as.data.frame(pivot_longer(sim_percent_long,
                                               cols = -Patient,
                                               names_to = "Cluster",
                                               values_to = "Percentage"))

comb <- rbind(do.call(rbind, cellgeo_percent),
              do.call(rbind, music_percent),
              do.call(rbind, DWLS_percent_NN),
              do.call(rbind, Lin_percent_NN))

comb$Data <- ifelse(grepl("Human", comb$Cluster), "Neuron", "Non_neuronal")

CCC_all_sub <- function(preds, obs, method, rep){
  sub <- preds[preds$Method == method & preds$Rep == rep, ]
  
  stopifnot(all(paste(sub$Patient, sub$Cluster) == 
                  paste(obs$Patient, obs$Cluster)))
  
  CCC(sub$Percentage, obs$Percentage)
}

sim_percent_long_neuron <- sim_percent_long[grepl("Human", sim_percent_long$Cluster), ]
comb_neuron <- comb[grepl("Human", comb$Cluster), ]

sim_percent_long_NN <- sim_percent_long[!grepl("Human", sim_percent_long$Cluster), ]
comb_NN <- comb[!grepl("Human", comb$Cluster), ]

CCC_comb_neuron <- list()

for(i in c("CellGeometry", "MuSiC")){
  for(x in paste("Rep", 1:3)){
    CCC_comb_neuron[[paste0(i, "_", gsub(" ", "", x))]] <- CCC_all_sub(comb_neuron, sim_percent_long_neuron, i, x)
  }
}

#saveRDS(CCC_comb_neuron, file = "CCC/brain_CCC_out_neuron.rds")

CCC_comb_res_neuron <- data.frame()

for(i in names(CCC_comb_neuron)){
  temp <- CCC_comb_neuron[[i]]$rho.c
  temp$Method <- gsub("_.*", "", i)
  temp$Rep <- gsub(".*_", "", i)
  temp$Rep <- gsub("Rep", "", temp$Rep)
  CCC_comb_res_neuron <- rbind(CCC_comb_res_neuron,
                               temp)
}

CCC_comb_res_neuron$Method <- factor(CCC_comb_res_neuron$Method,
                                     levels = c("CellGeometry",
                                                "MuSiC"))

CCC_comb_mean_neuron <- CCC_comb_res_neuron %>%
  group_by(Method) %>%
  summarise("Mean" = mean(est),
            "SD" = sd(est),
            "N" = n(),
            "SEM" = SD/sqrt(N))

CCC_comb_mean_neuron$Mean <- formatC(CCC_comb_mean_neuron$Mean, 3)
CCC_comb_mean_neuron$SEM <- formatC(CCC_comb_mean_neuron$SEM, 3)

# write.table(CCC_comb_mean_neuron[ , c("Method", "Mean", "SEM")],
#             file = "CCC/CCC_comb_res_neuron_mean.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

CCC_comb_NN <- list()

for(i in c("CellGeometry", "MuSiC", "DWLS", "LinDeconSeq")){
  for(x in paste("Rep", 1:3)){
    CCC_comb_NN[[paste0(i, "_", gsub(" ", "", x))]] <- CCC_all_sub(comb_NN, sim_percent_long_NN, i, x)
  }
}

#saveRDS(CCC_comb_NN, file = "CCC/brain_CCC_out_NN.rds")

CCC_comb_res_NN <- data.frame()

for(i in names(CCC_comb_NN)){
  temp <- CCC_comb_NN[[i]]$rho.c
  temp$Method <- gsub("_.*", "", i)
  temp$Rep <- gsub(".*_", "", i)
  temp$Rep <- gsub("Rep", "", temp$Rep)
  CCC_comb_res_NN <- rbind(CCC_comb_res_NN,
                           temp)
}

CCC_comb_res_NN$Method <- factor(CCC_comb_res_NN$Method,
                                 levels = c("CellGeometry",
                                            "MuSiC",
                                            "DWLS",
                                            "LinDeconSeq"))

CCC_comb_mean_NN <- CCC_comb_res_NN %>%
  group_by(Method) %>%
  summarise("Mean" = mean(est),
            "SD" = sd(est),
            "N" = n(),
            "SEM" = SD/sqrt(N))

CCC_comb_mean_NN$Mean <- formatC(CCC_comb_mean_NN$Mean, 3)
CCC_comb_mean_NN$SEM <- formatC(CCC_comb_mean_NN$SEM, 3)

# write.table(CCC_comb_mean_NN[ , c("Method", "Mean", "SEM")],
#             file = "CCC/CCC_comb_res_NN_mean.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

####summary table####

amp_mean <- read.delim("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/CCC/CCC_comb_res_mean.txt",
                       row.names = 1)

typist_mean <- read.delim("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Cell_typist/workstation2/Dirichlet/CCC/CCC_comb_res_mean.txt",
                          row.names = 1)

tabula_mean <- read.delim("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/CCC/CCC_comb_res_mean.txt",
                          row.names = 1)

amp_mean$Data <- "AMP"
typist_mean$Data <- "Cell Typist"
tabula_mean$Data <- "Tabula Sapiens"

CCC_comb_mean_neuron$Data <- "Human Brain Cell Atlas (Neuron)"
CCC_comb_mean_NN$Data <- "Human Brain Cell Atlas (Non-neuronal)"

mean_df <- rbind(amp_mean,
                 typist_mean,
                 tabula_mean,
                 CCC_comb_mean_neuron[ , c("Method", "Mean", "SEM", "Data")],
                 CCC_comb_mean_NN[ , c("Method", "Mean", "SEM", "Data")])

mean_df$CCC <- paste(mean_df$Mean, "±", mean_df$SEM)

# write.table(mean_df[ , c("Data", "Method", "CCC")], file = "CCC/mean_all.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

####Bland-Altman plots####

amp_CCC <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/AMP_CCC_out.rds")
typist_CCC <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Cell_typist/workstation2/Dirichlet/typist_CCC_out.rds")
tabula_CCC <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/Tabula_CCC_out.rds")

balt_all <- function(CCC_res, data){
  comb_balt <- data.frame()
  
  for(i in names(CCC_res)){
    temp <- CCC_res[[i]]$blalt
    temp$Method <- gsub("_.*", "", i)
    temp$Rep <- gsub(".*_", "", i)
    temp$Rep <- gsub("Rep", "", temp$Rep)
    comb_balt <- rbind(comb_balt,
                       temp)
  }
  
  comb_balt$Data <- data
  comb_balt
}

balt_all_df <- rbind(balt_all(amp_CCC, "AMP"),
                     balt_all(typist_CCC, "Cell Typist"),
                     balt_all(tabula_CCC, "Tabula Sapiens"),
                     balt_all(CCC_comb_neuron, "Human Brain Cell Atlas (Neuron)"),
                     balt_all(CCC_comb_NN, "Human Brain Cell Atlas (Non-neuronal)"))

balt_plot <- function(df, method, stat){
  sub <- df[df$Method == method, ]
  
  stat_sub <- stat[stat$Method == method, ]

  metric <- paste0("CCC = ",
                   stat$CCC)
  
  ggplot(sub,
         aes(x = mean, y = delta)) +
    rasterize(geom_point(alpha = 0.1), scale = 0.6, dpi = 500) +
    geom_hline(yintercept = mean(sub$delta),
               color = "blue") +
    geom_hline(yintercept = mean(sub$delta) + (1.96* sd(sub$delta)),
               color = "red", linetype = "dashed") +
    geom_hline(yintercept = mean(sub$delta) - (1.96* sd(sub$delta)),
               color = "red", linetype = "dashed") +
    labs(x = "Mean", y = "Pred - Obs",
         title = method, subtitle = metric) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"))
}

balt_plot_all <- function(data){
  
  sub <- balt_all_df[balt_all_df$Data == data, ]
  
  stat_sub <- mean_df[mean_df$Data == data, ]
  
  fig <- lapply(unique(sub$Method), function(x){
    balt_plot(sub[sub$Method == x, ], x, stat_sub[stat_sub$Method == x, ])
  })
  
  annotate_figure(ggarrange(plotlist = fig, ncol = length(unique(sub$Method))),
                  top = text_grob(label = data))
}

library(ggplot2)
library(ggrastr)
library(ggpubr)

pdf("CCC/AMP_balt.pdf", width = 12.5, height = 3.65)
balt_plot_all("AMP")
dev.off()

pdf("CCC/typist_balt.pdf", width = 12.5, height = 3.65)
balt_plot_all("Cell Typist")
dev.off()

pdf("CCC/tabula_balt.pdf", width = 9.375, height = 3.65)
balt_plot_all("Tabula Sapiens")
dev.off()

pdf("CCC/neuron_balt.pdf", width = 6.25, height = 3.65)
balt_plot_all("Human Brain Cell Atlas (Neuron)")
dev.off()

pdf("CCC/NN_balt.pdf", width = 12.5, height = 3.65)
balt_plot_all("Human Brain Cell Atlas (Non-neuronal)")
dev.off()





