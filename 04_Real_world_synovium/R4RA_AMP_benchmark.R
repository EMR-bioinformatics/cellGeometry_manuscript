# Deconvolution of R4RA synovium with AMP dataset

counts <- readRDS("R4RA_counts.rds")

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))

library(cellGeometry)
library(tictoc)

mk <- readRDS("cellmarkers_R4RA_rawcounts.rds")

#200 favourable from initial analysis - testing different ngroups
mk_update <- updateMarkers(mk, nsubclass = 200, bulkdata = counts)

cellgeo_out <- function(n, df){
  mk_update <- updateMarkers(mk_update, ngroup = n, bulkdata = counts)
  tic()
  out <- deconvolute(mk_update, df, count_space = TRUE,
                     convert_bulk = F,  weight_method = "equal")
  toc(log = TRUE)
  time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  out$Time <- time
  out
}

ngroup<- c(5, 10, 20, 25, 50, 100, 200)

cellgeo_orig_ngroup <- lapply(ngroup, function(x){
  cellgeo_out(x, counts)
})

names(cellgeo_orig_ngroup) <- paste("Ngroup =", ngroup)

# saveRDS(cellgeo_orig_ngroup, "Cellgeo_AMP_rawcounts_countT_eq200_ngroup.rds")

####MuSiC2####

library(SingleCellExperiment)
library(MuSiC, lib.loc =  "/home/rachel/R/x86_64-pc-linux-gnu-library/4.2")

sce <- SingleCellExperiment(list(counts = celseq_counts[ ,is.na(metadata$subclass) == FALSE]),
                            colData = metadata[is.na(metadata$subclass) == FALSE, ])


tic()
music2fit <- music_prop(bulk.mtx = counts,
                        sc.sce = sce,
                        clusters = "subclass",
                        samples = "sample")
toc(log = TRUE)
#331.919 sec elapsed
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
music2fit$time <- time

# saveRDS(music2fit, file = "music2_AMP_rawcounts.rds")

####DWLS####

load("Signature.rdata")

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
  if (inherits(solDWLS, 'try-error')) solDWLS <- rep(NA, ncol(S))
  allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
}
toc(log = TRUE) 
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

colnames(allCounts_DWLS) <- colnames(counts)
rownames(allCounts_DWLS) <- paste0("SC-", rownames(allCounts_DWLS))

table(is.na(allCounts_DWLS[1, ]))

# FALSE  TRUE 
# 36   147 

#save(time, allCounts_DWLS, file = "DWLS_AMP_rawcounts.rdata")

#none that exact zero output
fill <- colSums(is.na(allCounts_DWLS)) == nrow(allCounts_DWLS)
allCounts_DWLS[ , fill] <- 0

#save(time, allCounts_DWLS, file = "DWLS_AMP_rawcounts.rdata")

####LinDeconSeq####

#restart

library(LinDeconSeq)
library(tictoc)

markerRes <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/LinDeconSeq_signature.rds")

tic()
fractions <- deconSeq(counts, markerRes$sigMatrix$sig.mat)
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 

#save(time, fractions, file = "LinDeconSeq_AMPv2_rawcounts.rdata")




