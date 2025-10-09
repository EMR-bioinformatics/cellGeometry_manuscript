library(devtools)
setwd("/users/myles/documents/github/cellGeometry")
load_all()

library(plotly)
library(volcano3D)

load("/Users/myles/R/Deconv/AMP_scRNAseq_data_and_annotations.RData")  # single cell data
tpmdata <- readRDS("/Users/myles/R/R4RA/txi.R4RA.counts.tpm.final.rds")  # bulk RNA-Seq

# reformat metadata to avoid unclassified cells
metadata$type[metadata$type == "Empty"] <- NA
metadata$type <- factor(metadata$type,
                        levels = c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))
metadata$subclass <- factor(metadata$subclass, levels = col_order)

mk <- cellMarkers(celseq_counts,
                  bulkdata = tpmdata,
                  subclass = metadata$subclass,
                  cellgroup = metadata$type)

# rescale
df <- scaleSphere(mk$groupmeans[, c("Fibroblast", "Monocyte", "B cell")])
angle <- matrixStats::rowMins(acos(df)) * 180 / pi
df <- as.data.frame(df)
df$angle <- angle
df$max <- rowMaxs(mk$groupmeans[, c("Fibroblast", "Monocyte", "B cell")])
df <- df[complete.cases(df), ]
df <- df[order(df$max, decreasing = FALSE), ]
htext <- paste0(rownames(df), "<br>",
                "Expression ", formatC(df$max, digits = 3))
df$max[df$max > 8] <- 8
scheme <- hcl.colors(10, "Plasma", rev = TRUE)

# 3d sphere 
p <- plot_ly(df, x = ~Fibroblast, y = ~Monocyte, z = ~`B cell`,
        color = ~max, colors = scheme,
        mode = "markers", type = "scatter3d",
        marker = list(symbol = "circle", size = ~max,
                      sizeref = 0.02,
                      sizemin = 4, sizemode = "area", line = list(width = 0)),
        text = htext) %>%
  add_animation(speed = 720)

dd <- p$x$visdat[[1]]()

