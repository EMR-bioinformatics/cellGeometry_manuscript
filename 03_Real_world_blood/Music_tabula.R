####Background####

# Here using MuSiC to deconvolute PEAC blood using tabula blood

PEAC.blood <- readRDS("peac_bld_counts.rds")

library(zellkonverter)
library(SingleCellExperiment)

h5 <- readH5AD("/media/lvm1/tabula/4eb58518-23ce-43e8-89ca-9e9fa87c081d.h5ad",
               use_hdf5 = TRUE, reader = "R")

mat <- h5@assays@data$X
rownames(mat) <- rownames(h5)

dim(mat)
#61759 85233
#so cannot be used for DWLS and LinDeconSeq

library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

library(cellGeometry)

rownames(mat) <- gene2symbol(rownames(mat), ensDb_v110)

meta <- h5@colData@listData

library(MuSiC)

sce <- SingleCellExperiment(list(counts = mat),
                            colData = data.frame("subclass" = meta$cell_type,
                                                 "donor_id" = meta$donor_id))

library(tictoc)
tic()
music2fit <- music_prop(bulk.mtx = PEAC.blood,
                        sc.sce = sce,
                        clusters = "subclass",
                        samples = "donor_id")
toc(log = TRUE)
#
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
music2fit$time <- time
#saveRDS(music2fit, file = "music2_PEAC_tabula.rds")


