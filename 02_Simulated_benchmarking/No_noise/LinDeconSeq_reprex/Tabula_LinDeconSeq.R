# Generating Tabula Sapiens signature using LinDeconSeq for deconvolution

library(zellkonverter)
library(SingleCellExperiment)

h5 <- readH5AD("/media/lvm1/tabula/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad",
               use_hdf5 = TRUE, reader = "R")
#scRNA-seq dataset available from:  https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5

mat <- h5@assays@data$X
rownames(mat) <- rownames(h5)

meta <- h5@colData@listData

slim <- function(subclass, maxlen = 200) {
  tab <- table(subclass)
  ok <- tab >= 10 & tab <= maxlen
  sam <- which(tab > maxlen)
  ind <- which(subclass %in% levels(subclass)[ok])
  if (length(sam) > 0) {
    ind2 <- unlist(lapply(sam, function(i) {
      subcl <- which(subclass == levels(subclass)[i])
      sample(subcl, maxlen)
    }))
    ind <- sort(c(ind, ind2))
  }
  ind
}

ind <- slim(meta$cell_type, 100)
# generating index to decrease the size of the scRNA-seq data
# so that matrix is small enough to be read by LinDeconSeq

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

meta_sub <- data.frame("Cell_type" = meta$cell_type[ind],
                       "Cell" = meta$observation_joinid[ind])

length(unique(meta_sub$Cell)) == nrow(meta_sub)
#TRUE

meta_sub$Cell <- paste0("Cell", 1:nrow(meta_sub))

library(tidyverse)
phes <- as.data.frame(meta_sub %>% pivot_wider(Cell_type,
                                               names_from = Cell,
                                               values_from = Cell))
rownames(phes) <- phes$Cell_type
phes <- phes[ , -1]
phes <- as.matrix(phes)
phes <- ifelse(is.na(phes), 0, 1)

#table(colSums(phes))
#all 1s so fine

mat_sub <- as.matrix(mat[ ,ind])

colnames(mat_sub) <- paste0("Cell", 1:nrow(meta_sub))

library(LinDeconSeq)

markerRes <- findMarkers(mat_sub,
                         phes = phes)

# [INFO] 15739 samples and 61759 genes in the reference profile
# [INFO] 1836 candidate cell type-specific genes detected with q.cut 0.01
# [INFO] Select seed genes and random permutating...
# [INFO] Assigning cell type-specific genes to each cell subset
# [INFO] 0 low confidence marker genes have been filtered out.                                                                                                         
# [INFO] Optimizing cell type-specific genes to derive signature matrix...
# Error in if (sum(tmp.lab) > 1) { : missing value where TRUE/FALSE needed




