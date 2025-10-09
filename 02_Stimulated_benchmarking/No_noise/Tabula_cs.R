####Background####

# Tabula Sapiens heatmap of the cosine similiarity 

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2")

mk <- readRDS("tabula_markers_dualmeans_rerun.rds")

library(zellkonverter)
library(SingleCellExperiment)

h5 <- readH5AD("/media/lvm1/tabula/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad",
               use_hdf5 = TRUE, reader = "R")

meta <- h5@colData@listData

cellgroup <- data.frame(subclass = meta$cell_type, 
                        group = meta$broad_cell_class)

cellgroup <- distinct(cellgroup)

cellgroup$group <- as.vector(cellgroup$group)

celltable <- as.data.frame(mk$cell_table)
colnames(celltable) <- "Group"

celltable$Fine_group <- cellgroup$group[match(rownames(celltable),
                                              cellgroup$subclass)]

celltable$Groupv2 <- as.vector(celltable$Group)
celltable$Groupv2[celltable$Group == "Immune"] <- celltable$Fine_group[celltable$Group == "Immune"]

library(stringr)

celltable$Groupv2 <- str_to_sentence(celltable$Groupv2)
celltable$Groupv2 <- gsub(" cell", "", celltable$Groupv2)
celltable$Groupv2[celltable$Groupv2 == "Lymphocyte of b lineage"] <- "B cells"
celltable$Groupv2[celltable$Groupv2 == "Stem"] <- "Stem cells"
celltable$Groupv2[celltable$Groupv2 == "T"] <- "T cells"
celltable$Groupv2[celltable$Groupv2 == "Erythroid lineage"] <- "Erythroid"

rownames(celltable)[rownames(celltable) == "CD4-positive, alpha-beta T cell"] <- "CD4 T cell"
rownames(celltable)[rownames(celltable) == "naive thymus-derived CD4-positive, alpha-beta T cell"] <- "naive thymus-derived CD4 T cell"

colscheme_orig <- c("Germline" = "black",
                    "Stromal" = "green3",
                    "Endothelium" = "darkorange",
                    "Epithelium" = "thistle",
                    "Neural" = "bisque",
                    "Dendritic" = "sienna4",
                    "Granulocyte" = "gold",
                    "Erythroid" = "indianred",
                    "Hematopoietic" = "brown1",
                    "Myeloid leukocyte" = "firebrick3",
                    "Innate lymphoid" = "aquamarine",
                    "B cells" = "steelblue",
                    "T cells" = "purple2",
                    "Stem cells" = "grey60")

colscheme <- c("#f34e5a",
               "#59ec8b",
               "#d74eb5",
               "#008b6a",
               "#b0004e",
               "#00abf7",
               "#dd6f22",
               "#0060c4",
               "#f4d06b",
               "#a17df8",
               "#967100",
               "#535c9f",
               "#7d0044",
               "#d49ee0")

names(colscheme) <- names(colscheme_orig)

highlight <- c("CD4 T cell",
               "naive thymus-derived CD4 T cell",
               "fibroblast",
               "monocyte",
               "adventitial cell",
               "thymic fibroblast type 2",
               "intermediate monocyte",
               "classical monocyte")

cs_plot <- function(nsub){
  mk <- updateMarkers(mk, nsubclass = nsub)
  
  tabula_cs <- cos_similarity(mk)
  
  colnames(tabula_cs)[colnames(tabula_cs) == "CD4-positive, alpha-beta T cell"] <- "CD4 T cell"
  rownames(tabula_cs)[rownames(tabula_cs) == "CD4-positive, alpha-beta T cell"] <- "CD4 T cell"
  colnames(tabula_cs)[colnames(tabula_cs) == "naive thymus-derived CD4-positive, alpha-beta T cell"] <- "naive thymus-derived CD4 T cell"
  rownames(tabula_cs)[rownames(tabula_cs) == "naive thymus-derived CD4-positive, alpha-beta T cell"] <- "naive thymus-derived CD4 T cell"
 
  celltable <- celltable[match(rownames(tabula_cs),
                               rownames(celltable)), ]
  
  celltable$Groupv2 <- factor(celltable$Groupv2,
                              levels = names(colscheme))
  
  stopifnot(rownames(celltable) == rownames(tabula_cs))
  
  annotate_label <- rowAnnotation("Cell_group" = celltable$Groupv2,
                                  col = list(Cell_group = c(colscheme)),
                                  foo = anno_mark(at = which(rownames(tabula_cs) %in% highlight),
                                                  labels = rownames(tabula_cs)[which(rownames(tabula_cs) %in% highlight)],
                                                  labels_gp = gpar(fontsize = 9)),
                                  show_annotation_name = FALSE,
                                  annotation_legend_param = list(
                                    Cell_group = list(title = 'cell group',
                                                      title_gp = gpar(fontface = 1,
                                                                      fontsize = 6),
                                                      labels_gp = gpar(fontsize = 6))))
  ht1 <- Heatmap(tabula_cs, 
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 right_annotation = annotate_label,
                 show_row_dend = FALSE,
                 heatmap_legend_param = list(title = 'cosine similarity',
                                             title_gp = gpar(fontface = 1,
                                                             fontsize = 6),
                                             labels_gp = gpar(fontsize = 6)))
  
  draw(ht1, heatmap_legend_side = "right")
}

pdf("cs_heatmap/tabula_cosine_50v2.pdf", width = 7.5, height = 4.6)
cs_plot(50)
dev.off()






