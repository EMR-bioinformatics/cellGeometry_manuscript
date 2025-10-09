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
cellgrp <- meta$Majority_voting_CellTypist_high

# reduce dataset to only blood (optional)
subcl[meta$tissue != "blood"] <- NA
cellgrp[meta$tissue != "blood"] <- NA

mk <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp,
                  dual_mean = TRUE, cores = 8)


library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]
mk <- gene2symbol(mk, ensDb_v110)

# simulated bulk
set.seed(3)
sim_counts <- generate_samples(mk, 25)
sim_percent <- sim_counts / rowSums(sim_counts) * 100
sim_pseudo <- simulate_bulk(mk, sim_counts)

sim_sampled <- simulate_bulk(mat, sim_counts, subcl, times = 3)
rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)

fit2 <- deconvolute(mk, sim_sampled,
                    use_filter = FALSE, arith_mean = TRUE)

mkv2 <- updateMarkers(mk,
                      nsubclass = 5)

#for tidy heatmap shorten cell group names

temp <- as.vector(mkv2[["cell_table"]])
names(temp) <- names(mkv2[["cell_table"]])
temp[temp == "Monocytes"] <- "Mono"
temp[temp == "Cycling cells"] <- "Cycling"
temp[temp == "Early MK"] <- "MK"
temp[temp == "Macrophages"] <- "Mph"
temp[temp == "Mast cells"] <- "Mast"
temp[temp == "Plasma cells"] <- "PC"
temp <- factor(temp,
               levels = c("ILC", "Mono", "Cycling", "T cells",
                          "DC", "MK", "HSC/MPP", "Mph",
                          "Mast", "B cells", "PC", "pDC"))

mkv2[["cell_table"]] <- temp

# gene signature
pdf(file = "typist_sig.pdf", height = 10, width = 6)
signature_heatmap(mkv2)
dev.off()

# spillover
pdf(file = "spillover_notext.pdf", height = 5.25, width = 5)
spillover_heatmap(mkv2, text = FALSE)
dev.off()

# compensation
comp_heatmap(fit2, text = FALSE)
#550x550 

