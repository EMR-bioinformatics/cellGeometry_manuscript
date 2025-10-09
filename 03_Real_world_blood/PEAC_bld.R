# testing deconvolution of PEAC cell report blood with cell typist

counts <- readRDS("peac_bld_counts.rds")

####MuSiC2####

typist <- readRDS("/media/lvm1/celltypist/2ac906a5-9725-4258-8e36-21a9f6c0302a.rds")
meta <- typist@meta.data

subcl <- meta$Majority_voting_CellTypist
subcl[meta$tissue != "blood"] <- NA
cellgrp <- meta$Majority_voting_CellTypist_high
cellgrp[meta$tissue != "blood"] <- NA

library(SingleCellExperiment)
library(MuSiC)
library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

mat <- typist@assays$RNA$counts

rownames(mat) <- gene2symbol(rownames(mat), ensDb_v110)

sce <- SingleCellExperiment(list(counts = mat[ , is.na(subcl) == FALSE]),
                            colData = data.frame("subclass" = subcl[is.na(subcl) == FALSE],
                                                 "donor_id" = meta$donor_id[is.na(subcl) == FALSE]))


tic()
music2fit <- music_prop(bulk.mtx = counts,
                        sc.sce = sce,
                        clusters = "subclass",
                        samples = "donor_id")
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
music2fit$time <- time

# saveRDS(music2fit, file = "music2_PEAC_typist_rawcounts.rds")

####DWLS####

library(DWLS)

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Cell_typist/workstation2/DWLS_signature.rdata")

allCounts_DWLS<-NULL

tic() 
for(j in 1:(dim(counts)[2])){
  S <- Signature
  Bulk<-counts[, j]
  names(Bulk)<- rownames(counts)
  Genes <- intersect(rownames(S), names(Bulk))
  B <- Bulk[Genes]
  S <- S[Genes,]
  solDWLS<-try(solveDampenedWLS(S,B), silent = TRUE)
  if (inherits(solDWLS, 'try-error')) solDWLS <- rep(NA, ncol(S)) #NA added here
  allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
}
toc(log = TRUE) 
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

colnames(allCounts_DWLS) <- colnames(counts)

# save(time, allCounts_DWLS, file = "DWLS_PEAC_typist_rawcounts.rdata")

####LinDeconSeq####

# restart

library(LinDeconSeq)

Signature <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Cell_typist/workstation2/LinDeconSeq_signature.rds")

library(tictoc)
tic()
fractions <- deconSeq(counts, Signature$sigMatrix$sig.mat)
toc(log = TRUE) 
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

# save(time, fractions, file = "LinDeconSeq_PEAC_typist.rdata")


