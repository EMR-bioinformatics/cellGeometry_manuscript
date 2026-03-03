# Downsizing tabula dataset for DWLS deconvolution 

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2")

library(zellkonverter)
library(SingleCellExperiment)

h5 <- readH5AD("/media/lvm1/tabula/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad",
               use_hdf5 = TRUE, reader = "R")

mat <- h5@assays@data$X
rownames(mat) <- rownames(h5)

meta <- h5@colData@listData

# slim <- function(subclass, maxlen = 200) {
#   tab <- table(subclass)
#   ok <- tab >= 10 & tab <= maxlen
#   sam <- which(tab > maxlen)
#   ind <- which(subclass %in% levels(subclass)[ok])
#   if (length(sam) > 0) {
#     ind2 <- unlist(lapply(sam, function(i) {
#       subcl <- which(subclass == levels(subclass)[i])
#       sample(subcl, maxlen)
#     }))
#     ind <- sort(c(ind, ind2))
#   }
#   ind
# }
# 
# ind <- slim(meta$cell_type, 100) #200 R was grumbling downstream

ind <- readRDS("tabula_subset_index.rds")

celltypedf <- data.frame("Orig" = unique(as.vector(meta$cell_type)))

library(checkmate)
check_names(unique(as.vector(meta$cell_type)), "strict")
#"Must have names according to R's variable naming conventions, but element 1 does not comply"

meta$cell_type <- gsub("-", "_", meta$cell_type)

check_names(unique(as.vector(meta$cell_type)), "strict")
#"Must have names according to R's variable naming conventions, but element 1 does not comply"

meta$cell_type <- gsub(",", " ", meta$cell_type)

check_names(unique(as.vector(meta$cell_type)), "strict")
#"Must have names according to R's variable naming conventions, but element 1 does not comply"

meta$cell_type <- gsub("  ", " ", meta$cell_type)

check_names(unique(as.vector(meta$cell_type)), "strict")
#"Must have names according to R's variable naming conventions, but element 1 does not comply"

meta$cell_type <- gsub(" ", "_", meta$cell_type)
check_names(unique(as.vector(meta$cell_type)), "strict")
#"Must have names according to R's variable naming conventions, but element 132 does not comply"

meta$cell_type <- gsub("\\+", "plus", meta$cell_type)
check_names(unique(as.vector(meta$cell_type)), "strict")
#TRUE

celltypedf$convert <- unique(as.vector(meta$cell_type))

#saveRDS(celltypedf, file = "DWLS_celltype_convert.rds")

subcl <- meta$cell_type[ind]
mat_sub <- as.matrix(mat[ ,ind])

rm(h5)
rm(mat)
rm(meta)

colnames(mat_sub) <- paste0("Cell", 1:length(subcl))

subcl <- factor(subcl)

#dim(mat_sub)
#61759 15739

rm <- which(apply(mat_sub, 1, function(x) all(x==0)))
#2628

mat_sub <- mat_sub[!rownames(mat_sub) %in% names(rm), ] # remove genes always zero
# dim(mat_sub)
# 59131 15739

rm2 <- which(apply(mat_sub, 1, function(x) length(which(x>0))<10))
#12910

mat_sub <- mat_sub[!rownames(mat_sub) %in% names(rm2), ]
# remove genes where only 10 out of 15739 cells express these genes
#dim(mat_sub)
#46221 15739

# save(mat_sub, subcl, 
#      file = "Tabula_subset.rdata")