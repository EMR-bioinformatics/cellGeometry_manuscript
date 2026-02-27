# Generating Brain neurons signature using LinDeconSeq for deconvolution

library(zellkonverter)
library(SingleCellExperiment)

brain <- readH5AD("/media/lvm1/brain/c2f66cd5-4ff4-4578-876c-55783a57cf8f.h5ad",
                  use_hdf5 = TRUE, reader = "R")
#Neuron scRNA-seq dataset available from: https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443

mat <- brain@assays@data$X
rownames(mat) <- rownames(brain)  
meta <- brain@colData@listData

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

ind <- slim(meta$roi, 100) 

library(checkmate)
check_names(unique(as.vector(meta$roi)), "strict")
#[1] "Must have names according to R's variable naming conventions, but element 1 does not comply"

meta$roi <- gsub(" ", "_", meta$roi)
check_names(unique(as.vector(meta$roi)), "strict")
#[1] "Must have names according to R's variable naming conventions, but element 4 does not comply"

meta$roi <- gsub("-", "_", meta$roi)
check_names(unique(as.vector(meta$roi)), "strict")
#TRUE

meta_sub <- data.frame("Cell_type" = meta$roi[ind],
                       "Cell" = meta$observation_joinid[ind])

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

# [INFO] 10444 samples and 59236 genes in the reference profile
# [INFO] 2899 candidate cell type-specific genes detected with q.cut 0.01
# [INFO] Select seed genes and random permutating...
# [INFO] Assigning cell type-specific genes to each cell subset
# [INFO] 0 low confidence marker genes have been filtered out.                                                                                                         
# [INFO] Optimizing cell type-specific genes to derive signature matrix...
# Error in if (sum(tmp.lab) > 1) { : missing value where TRUE/FALSE needed




