# deconvoluting STRAP synovium (AMP dataset)

load("/media/gcpeac/STRAP/DB_releases/v202202/STRAP_db_v202202.RData")

mk <- readRDS("cellmarkers.rds") #no filtering of gene names 

rnasamples <- STRAP_db_v202202$RNAseq_samples
rnasamples <- rnasamples[rnasamples$Outliers == FALSE, ]
rnasamplesBL <- rnasamples[rnasamples$Visit == 2, ]
strap_meta <- STRAP_db_v202202$Outcome

strap_metaBl <- strap_meta[which((strap_meta$Patient.I.D.%in%rnasamplesBL $anonymised_PatientID)&
                                   (strap_meta$Visit==3)),]
remove <- setdiff(rnasamples$anonymised_PatientID[rnasamples$Visit == 2], strap_metaBl$Patient.I.D.)
#2 with no baseline clinical data

rnasamplesv2 <- rnasamples[-which(rnasamples$anonymised_PatientID %in% remove &
                                    rnasamples$Visit == 2), ]

counts <- STRAP_db_v202202$RNAseq_data$txi$counts
counts <- counts[ , rnasamplesv2$anonymised_Sample_Name]

mk_update <- updateMarkers(mk, nsubclass = 200, ngroup = 5, bulkdata = counts)

tic()
cellgeo <- deconvolute(mk_update, counts, count_space = TRUE,
                       convert_bulk = F, weight_method = "equal",
                       comp_amount = 1)
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
cellgeo$Time <- time

#saveRDS(cellgeo, "Cellgeo_AMP_rawcounts_countT_eq200.rds")

mk_updatev2 <- updateMarkers(mk, nsubclass = 200, ngroup = 10, bulkdata = counts)

tic()
cellgeov2 <- deconvolute(mk_updatev2, counts, count_space = TRUE,
                         convert_bulk = F, weight_method = "equal",
                         comp_amount = 1)
toc(log = TRUE)
time <- as.numeric(gsub(" sec.*", "", unlist(tic.log(format = TRUE))))
tic.clearlog() 
cellgeov2$Time <- time

#saveRDS(cellgeov2, "Cellgeo_AMP_rawcounts_countT_eq200_ngroup10.rds")

