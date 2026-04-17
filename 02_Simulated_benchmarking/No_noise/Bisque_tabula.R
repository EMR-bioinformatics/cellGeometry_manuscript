####Tabula####

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/Bisque_tabula")
library(zellkonverter)
library(SingleCellExperiment)
library(Biobase)
library(tictoc)
library(BisqueRNA)
library(AnnotationHub)
library(cellGeometry)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

#needs to be in Expression Set format
#cannot accept delayed matrix 
load("../Tabula_subset.rdata") #same as DWLS and LinDeconSeq
load("../tabula_simulated_workstation2.rdata")

rownames(mat_sub) <- gene2symbol(rownames(mat_sub), ensDb_v110)

sc.pheno <- data.frame(check.names = F, check.rows = F,
                       stringsAsFactors = F,
                       row.names = colnames(mat_sub),
                       SubjectName = donor_id,
                       subclass = subcl)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "subclass"),
                      row.names=c("SubjectName",
                                  "subclass"))

sc.pdata <- new("AnnotatedDataFrame",
                data = sc.pheno,
                varMetadata = sc.meta)

sc.eset <- Biobase::ExpressionSet(assayData = mat_sub,
                                  phenoData = sc.pdata)

bisque_out <- function(times, rep){
  samples_sim <- sim_sampled[[paste("Times", times)]][[paste("Rep", rep)]]
  
  samples_sim <- ExpressionSet(assayData = samples_sim)
  
  tic()
  out <- ReferenceBasedDecomposition(bulk.eset = samples_sim, 
                                     sc.eset = sc.eset,
                                     use.overlap = FALSE,
                                     cell.types = "subclass")
  toc(log = TRUE)
  mk_time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
  tic.clearlog() 
  
  list(output = out,
       time = mk_time)
}

bisque <- list()

for(x in paste("Rep", 1:5)){
  bisque[["Times 3"]][[x]] <- bisque_out(times = 3,
                                         rep = as.numeric(gsub("Rep ", "", x)))
}

saveRDS(bisque, file = "bisque_tabula_output.rds")