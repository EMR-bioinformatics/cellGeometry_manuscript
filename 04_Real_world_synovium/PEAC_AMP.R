# deconvoluting PEAC synovium (AMP dataset)

load("/media/gcpeac/Rachel/Merging_datasets/RNA_BL_2nd/merged4David.RData")

PEAC_meta <- meta[meta$Cohort == "PEAC1", ]

PEAC_counts <- rawcounts[ , colnames(rawcounts) %in% PEAC_meta$SeqID]

library(cellGeometry)
library(tictoc)

mk <- readRDS("cellmarkers.rds")

tic()
cellgeo <- deconvolute(mk_update, PEAC_counts, count_space = TRUE,
                       convert_bulk = F, weight_method = "equal",
                       comp_amount = 1)
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
cellgeo$Time <- time

#saveRDS(cellgeo, "PEAC1_AMP_rawcounts_countT_eq200.rds")

mk_updatev2 <- updateMarkers(mk, nsubclass = 200, ngroup = 10, bulkdata = PEAC_counts)

tic()
cellgeov2 <- deconvolute(mk_updatev2, PEAC_counts, count_space = TRUE,
                         convert_bulk = F, weight_method = "equal",
                         comp_amount = 1)
toc(log = TRUE)
#5.239 sec elapsed
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
cellgeov2$Time <- time

#saveRDS(cellgeov2, "PEAC1_AMP_rawcounts_countT_eq200_ngroup10.rds")
