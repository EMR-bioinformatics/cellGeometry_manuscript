####Background####

# Generating blood signature for cellGeometry and then deconvoluting

# initial regroup for cell typist 

library(zellkonverter)
library(SingleCellExperiment)
library(cellGeometry)

typist_h5 <- readH5AD("/media/lvm1/celltypist/2ac906a5-9725-4258-8e36-21a9f6c0302a.h5ad",
                      use_hdf5 = TRUE, reader = "R")

mat <- typist_h5@assays@data$X
rownames(mat) <- rownames(typist_h5)
meta <- typist_h5@colData@listData

table(meta$Majority_voting_CellTypist)

subcl <- meta$Majority_voting_CellTypist
cellgrp <- as.vector(meta$Majority_voting_CellTypist_high)

# reduce dataset to only blood
subcl[meta$tissue != "blood"] <- NA
cellgrp[meta$tissue != "blood"] <- NA

temp <- table(cellgrp)
cellgrp[cellgrp == "Cycling cells"] <- "T cells"
cellgrp[cellgrp %in% c("Early MK", "HSC/MPP", "Promyelocytes")] <- "Progenitor/misc"
cellgrp[cellgrp %in% c("Macrophages", "pDC", "DC")] <- "Monocytes" 
cellgrp[cellgrp == "Plasma cells"] <- "B cells" 
cellgrp[cellgrp == "Mast cells"] <- "granulocyte" 

table(cellgrp)

library(cellGeometry)
library(tictoc)

tic()
celltypist <- cellMarkers(mat, subclass = subcl,
                          cellgroup = cellgrp,
                          remove_subclass = c("Helper T cells", "Cytotoxic T cells"),
                          dual_mean = T, cores = 8)

toc(log = TRUE)
mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog()

celltypist$mk_time <- mk_time

library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

celltypist  <- gene2symbol(celltypist, ensDb_v110)

# saveRDS(celltypist, "celltypist_regroup_hd5.rds")

signature_heatmap(celltypist, top = 5)

# retrieve tabula (blood)

h5 <- readH5AD("/media/lvm1/tabula/4eb58518-23ce-43e8-89ca-9e9fa87c081d.h5ad",
               use_hdf5 = TRUE, reader = "R")

mat <- h5@assays@data$X
rownames(mat) <- rownames(h5)

meta <- h5@colData@listData

tic()
tabula <- cellMarkers(mat, subclass = meta$cell_type, 
                  cellgroup = meta$broad_cell_class, 
                  dual_mean = TRUE, cores = 8)

toc(log = TRUE)
tabula_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

tabula$tabula_time <- tabula_time
tabula <- gene2symbol(tabula, ensDb_v110)

# saveRDS(tabula, file = "tabula_bld_markers_dualmeans.rds")

signature_heatmap(tabula, top = 5)

# merge

gp <- c(names(tabula$subclass_table))
except <- setdiff  # use with pipe

filter <- gp|> except(c("basophil", "neutrophil"))

merge_mk_hd5 <- mergeMarkers(celltypist, tabula,
                             remove_subclass = filter,
                             transform = "qq")

#saveRDS(merge_mk_hd5, "merge_markers_hd5.rds")

signature_heatmap(merge_mk_hd5, top = 5)

merge_mk_hd5_consol <- collapse_group(merge_mk_hd5, groups = c("granulocyte", "granulocyte.1"))

#saveRDS(merge_mk_hd5_consol, "merge_markers_hd5_collapse.rds")

signature_heatmap(merge_mk_hd5_consol, top = 5)

####deconvolution####

PEAC.blood <- readRDS("peac_bld_counts.rds")

mk_update <- updateMarkers(merge_mk_hd5_consol,
                           #ngroup = 5,
                           bulk = PEAC.blood)

tic()
cellgeo <- deconvolute(mk_update, countsv2)
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog()
cellgeo$Time <- time

#saveRDS(cellgeo, "Cellgeo_PEAC_merge.rds")














