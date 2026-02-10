####Background####

# This script is a reprex in deconvoluting R4RA synovial samples using 
# DWLS and the AMP1 database 

####Signature curation####

load("AMP_scRNAseq_data_and_annotations.RData")
# AMP1 synovium scRNA-seq 
# Available from immPort accession code SDY998
# https://www.immport.org/shared/study/SDY998/summary 

celseq_counts <- as.matrix(celseq_counts)

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

metadata$subclass <- gsub("SC-", "", metadata$subclass)
col_order <- c(paste0('F', 1:4), paste0('M', 1:4), paste0('T', 1:6),
               paste0('B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

dataSC <- celseq_counts[ , is.na(metadata$subclass) == FALSE]
labels <- metadata$subclass[is.na(metadata$subclass) == FALSE]

library(DWLS)
# CRAN installed https://github.com/sistia01/DWLS
# Pipeline/functions guide https://github.com/dtsoucas/DWLS

Signature <- buildSignatureMatrixMAST(scdata = dataSC,
                                      id = labels,
                                      path = "results")

#save(Signature, time, file = "Signature.rdata")

# to save time for users, Signature.rdata is available with the reprex

####Initial attempt of deconvolution####

counts <- readRDS("R4RA_counts_anon.rds")
# R4RA bulk RNA-seq gene expression matrix 

allCounts_DWLS <- NULL

for(j in 1:(dim(counts)[2])){
  S <- Signature
  Bulk <- counts[, j]
  names(Bulk) <- rownames(counts)
  Genes <- intersect(rownames(S), names(Bulk))
  B <- Bulk[Genes]
  S <- S[Genes,]
  solDWLS <- solveDampenedWLS(S,B)
  allCounts_DWLS <- cbind(allCounts_DWLS,solDWLS)
}

# Error in solve.QP(D, d, A, bzero) : 
#   constraints are inconsistent, no solution!

####Deconvolution (skip those fail)####

allCounts_DWLS <- NULL

for(j in 1:(dim(counts)[2])){
  S <- Signature
  Bulk <- counts[, j]
  names(Bulk) <- rownames(counts)
  Genes <- intersect(rownames(S), names(Bulk))
  B <- Bulk[Genes]
  S <- S[Genes,]
  solDWLS <- try(solveDampenedWLS(S,B), silent = TRUE)
  if (inherits(solDWLS, 'try-error')) solDWLS <- rep(NA, ncol(S))
  allCounts_DWLS <- cbind(allCounts_DWLS, solDWLS)
}

# saveRDS(allCounts_DWLS, "R4RA_AMP1_DWLS.rds")

table(colSums(is.na(allCounts_DWLS)) == nrow(allCounts_DWLS))
# FALSE  TRUE 
# 36   147 

# so 80.3% failure rate 








