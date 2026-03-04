####Background####

# benchmark plots of brain simulated data

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain")

library(cellGeometry)

cellgeo <- readRDS("cellgeo_dirichlet_output.rds")
music2 <- readRDS("music_dirichlet_output.rds")
music2_NN <- readRDS("music_dirichlet_output_NN.rds")
DWLS <- readRDS("DWLS_dirichlet_output_neuron.rds")
DWLS_NN <- readRDS("DWLS_dirichlet_output_NN.rds")
Lin_NN <- readRDS("Lin_dirichlet_output_NN.rds")

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/brain_simulated_merged.rdata")

arrange_metric <- function(list, method){
  df <- data.frame()
  for(i in names(list)){
    for(x in names(list[[i]])){
      temp <- as.data.frame(list[[i]][[x]]$metrics_percent)
      temp <- add_column(temp,
                         "Cluster" = rownames(temp),
                         .before = 1)
      temp$Times <- as.numeric(gsub("Times ", "", i))
      temp$Rep <- as.numeric(gsub("Rep ", "", x))
      df <- rbind(df,
                  temp)
    }
  }
  
  df$Method <- method
  df
}

cellgeo_metric <- arrange_metric(cellgeo, "CellGeometry")

#organise music data 

for(x in paste("Rep", 1:5)){
  propdf <- music2_NN[["Times 1"]][[x]]$output$Est.prop.weighted * 100
  sim_percentNN <- sim_percent_merge[ , !grepl("Human", colnames(sim_percent_merge))]
  
  common <- intersect(colnames(sim_percentNN), colnames(propdf))
  
  propdf <- propdf[ , colnames(propdf) %in% common]
  
  propdf <- as.data.frame(propdf)
  
  propdf <- propdf[ , match(colnames(sim_percentNN), colnames(propdf))]
  
  stopifnot(all(colnames(sim_percentNN) == colnames(propdf)))
  
  metrics_percent <- metric_set(sim_percentNN, propdf)
  
  music2_NN[["Times 1"]][[x]]$metrics_percent <- metrics_percent
  
}


for(x in paste("Rep", 1:3)){
  propdf <- music2[["Times 1"]][[x]]$output$Est.prop.weighted * 100
  sim_percent <- sim_percent_merge[ , grepl("Human", colnames(sim_percent_merge))]
  
  common <- intersect(colnames(sim_percent), colnames(propdf))
  
  propdf <- propdf[ , colnames(propdf) %in% common]
  
  missing <- setdiff(colnames(sim_percent), colnames(propdf))
  
  propdf <- as.data.frame(propdf)
  
  if(length(missing) > 0){
    for(i in missing){
      propdf[ , i] <- 0
    }
  }
  
  propdf <- propdf[ , match(colnames(sim_percent), colnames(propdf))]
  
  stopifnot(all(colnames(sim_percent) == colnames(propdf)))
  
  metrics_percent <- metric_set(sim_percent, propdf)
  
  music2[["Times 1"]][[x]]$metrics_percent <- metrics_percent
  
}

music2_metric <- rbind(arrange_metric(music2, "MuSiC"),
                       arrange_metric(music2_NN, "MuSiC"))

#arrange DWLS 
celltype <- readRDS("DWLS_neuron_convert.rds")

for(x in paste("Rep", 1:5)){
  propdf <- t(DWLS_NN[["Times 1"]][[x]]$output) * 100
  sim_percentNN <- sim_percent_merge[ , !grepl("Human", colnames(sim_percent_merge))]
  
  colnames(propdf) <- gsub("_", " ", colnames(propdf))
  
  propdf <- propdf[ , colnames(sim_percentNN)]
  
  stopifnot(all(colnames(sim_percentNN) == colnames(propdf)))
  
  metrics_percent <- metric_set(sim_percentNN, propdf)
  
  DWLS_NN[["Times 1"]][[x]]$metrics_percent <- metrics_percent
  
}

for(x in paste("Rep", 1:5)){
  propdf <- t(DWLS[["Times 1"]][[x]]$output) * 100
  sim_percent <- sim_percent_merge[ , grepl("Human", colnames(sim_percent_merge))]
  
  colnames(propdf) <- celltype$Orig[match(colnames(propdf),
                                          celltype$convert)]
  
  propdf <- propdf[ , colnames(sim_percent)]
  
  stopifnot(all(colnames(sim_percent) == colnames(propdf)))
  
  metrics_percent <- metric_set(sim_percent, propdf)
  
  DWLS[["Times 1"]][[x]]$metrics_percent <- metrics_percent
  
}

DWLS_metric <- rbind(arrange_metric(DWLS, "DWLS"),
                     arrange_metric(DWLS_NN, "DWLS"))


#arrange LinDeconSeq

for(x in paste("Rep", 1:5)){
  propdf <- Lin_NN[["Times 1"]][[x]]$output * 100
  sim_percentNN <- sim_percent_merge[ , !grepl("Human", colnames(sim_percent_merge))]
  
  colnames(propdf) <- gsub("_", " ", colnames(propdf))
  
  propdf <- propdf[ , colnames(sim_percentNN)]
  
  stopifnot(all(colnames(sim_percentNN) == colnames(propdf)))
  
  metrics_percent <- metric_set(sim_percentNN, propdf)
  
  Lin_NN[["Times 1"]][[x]]$metrics_percent <- metrics_percent
  
}

Lin_metric <- arrange_metric(Lin_NN, "LinDeconSeq")

metric_percent <- rbind(cellgeo_metric,
                        music2_metric,
                        DWLS_metric,
                        Lin_metric)

metric_percent <- metric_percent[metric_percent$Rep %in% c(1:3), ]

metric_percent$Method <- factor(metric_percent$Method,
                                levels = c("CellGeometry",
                                           "MuSiC",
                                           "DWLS",
                                           "LinDeconSeq"))

#saveRDS(metric_percent, "brain_metric_percent_all.rds")

write.table(metric_percent[metric_percent$Method != "LinDeconSeq", ],
            file = "/media/gcpeac/Rachel/Packages/cellGeometry_paper/Consolidate/Data_points/02_Simulated_benchmarking/brain_metrics.txt",
            append = FALSE, sep = "\t",
            col.names = NA, row.names = TRUE)

####non-neuronal metric####

mk <- readRDS("brain_cellmarkers_merged.rds")

cell_groups <- as.data.frame(mk$cell_table)
colnames(cell_groups) <- "Group"

metric_percent_NN <- metric_percent[!grepl("Human", metric_percent$Cluster), ]

metric_percent_NN$Celltype <- cell_groups$Group[match(metric_percent_NN$Cluster,
                                                      rownames(cell_groups))]

levels <-  c("Fibroblast", "Vascular", "Astrocyte", "Choroid plexus",
             "Oligodendrocyte precursor", "Oligodendrocyte",
             "Ependymal", "Bergmann glia",  "Microglia")

metric_percent_NN$Celltype <- factor(metric_percent_NN$Celltype,
                                     levels = levels)

library(ggplot2)

# independent y axis

Pv2 <- ggplot(metric_percent_NN, aes(x = Celltype, y = RMSE, 
                                     group = Celltype)) +
  geom_jitter(aes(color = Method), alpha = 0.8, width = 0.1)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "LinDeconSeq" = "#00B0F6",
                                "DWLS" = "#E76BF3"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 4, scale = "free_y")+
  labs(x = "", y = "RMSE") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))

P2v2 <- ggplot(metric_percent_NN, aes(x = Celltype, y = Rsq, 
                                      group = Celltype)) +
  geom_jitter(aes(color = Method), alpha = 0.8, width = 0.1)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "LinDeconSeq" = "#00B0F6",
                                "DWLS" = "#E76BF3"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 4, scale = "free_y")+
  labs(x = "", y = "Rsq") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))

ggarrange(Pv2,P2v2, ncol = 1, nrow = 2,
          align = "hv")
#850x700

# shared y axis

# RMSE
library(ggbreak)

ggplot(data = metric_percent_NN) +
  geom_jitter(aes(x = Celltype, y = RMSE, group = Celltype, color = Method),
              alpha = 0.8, 
              width = 0.2)+
  geom_boxplot(aes(x = Celltype, y = RMSE), 
               alpha = 0.5,
               outlier.shape = NA) + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "LinDeconSeq" = "#00B0F6",
                                "DWLS" = "#E76BF3"),
                     guide = "none") +
  labs(x = "", y = "RMSE") +
  scale_y_break(c(40, 40), ticklabels = c(40, 60), scale = 0.2) +
  facet_wrap(~Method, ncol = 4)+
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        axis.title.y = element_text(hjust = 0.65),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))
#800x350

####neuron metric####

metric_percent_neuron <- metric_percent[grepl("Human", metric_percent$Cluster), ]

metric_percent_neuron$Celltype <- as.vector(cell_groups$Group[match(metric_percent_neuron$Cluster,
                                                                    rownames(cell_groups))])

#Human A35r is miscellaneous
levels2 <- c(sort(unique(metric_percent_neuron$Celltype[metric_percent_neuron$Celltype != "Miscellaneous"])))

metric_percent_neuron$Celltype <- factor(metric_percent_neuron$Celltype,
                                         levels = levels2)

metric_percent_neuron <- metric_percent_neuron[is.na(metric_percent_neuron$Celltype) == FALSE, ]

# independent y axis

P3v2 <- ggplot(metric_percent_neuron, aes(x = Celltype, y = RMSE, 
                                          group = Celltype)) +
  geom_jitter(aes(color = Method), alpha = 0.6, width = 0.1)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "DWLS" = "#E76BF3"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 3, scale = "free_y")+
  labs(x = "", y = "RMSE") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        plot.margin = margin(0.2, 0.2, 0.2, 1.2, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

P4v2 <- ggplot(metric_percent_neuron, aes(x = Celltype, y = Rsq, 
                                          group = Celltype)) +
  geom_jitter(aes(color = Method), alpha = 0.6, width = 0.1)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "DWLS" = "#E76BF3"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 3, scale = "free_y")+
  labs(x = "", y = "Rsq") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        plot.margin = margin(0.2, 0.2, 0.2, 1.2, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

ggarrange(P3v2,P4v2, ncol = 1, nrow = 2,
          align = "hv")
#850x700

# shared y axis

# RMSE

ggplot(data = metric_percent_neuron) +
  geom_jitter(aes(x = Celltype, y = RMSE, group = Celltype, color = Method),
              alpha = 0.8, 
              width = 0.2)+
  geom_boxplot(aes(x = Celltype, y = RMSE), 
               alpha = 0.5,
               outlier.shape = NA) + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "DWLS" = "#E76BF3"),
                     guide = "none") +
  labs(x = "", y = "RMSE") +
  facet_wrap(~Method, ncol = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        plot.margin = margin(0.2, 0.2, 0.2, 1.4, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))
#800x350

####scatter plots####

######neuron#####
sim_percent <- sim_percent_merge[ , grepl("Human", colnames(sim_percent_merge))]
cellgeo_percent <- as.data.frame(cellgeo$`Times 1`$`Rep 1`$output$subclass$percent)

cellgeo_percent_neuron <- cellgeo_percent[ , grepl("Human", colnames(cellgeo_percent))]

#all(colnames(sim_percent) == colnames(cellgeo_percent_neuron))
#TRUE

plot_pred(sim_percent, as.matrix(cellgeo_percent_neuron), mk) + labs(x = "Observed (%)", y = "Predicted (%)")
#500x500

prop_arrange <- function(df){
  common <- intersect(colnames(sim_percent), colnames(df))
  propdf <- df[ , colnames(df) %in% common]
  
  missing <- setdiff(colnames(sim_percent), colnames(propdf))
  
  propdf <- as.data.frame(propdf)
  
  if(length(missing) > 0){
    for(i in missing){
      propdf[ , i] <- 0
    }
  }
  
  propdf <- propdf[ , match(colnames(sim_percent), colnames(propdf))]
  propdf
}

music_percent <- prop_arrange(music2$`Times 1`$`Rep 1`$output$Est.prop.weighted * 100)

all(colnames(sim_percent) == colnames(music_percent))
#TRUE

plot_pred(sim_percent, as.matrix(music_percent), mk, ellipse = 2) + labs(x = "Observed (%)", y = "Predicted (%)")

DWLS_percent <- t(DWLS$`Times 1`$`Rep 1`$output) * 100
colnames(DWLS_percent) <- celltype$Orig[match(colnames(DWLS_percent),
                                              celltype$convert)]
DWLS_percent <- DWLS_percent[ , colnames(sim_percent)]

all(colnames(sim_percent) == colnames(DWLS_percent))
#TRUE

plot_pred(sim_percent, as.matrix(DWLS_percent), mk, ellipse = 2) + labs(x = "Observed (%)", y = "Predicted (%)")
#550x550

####obtain signatures for heatmap####

ht <- signature_heatmap(mk, scale = "sphere", top = 5, text = FALSE, use_raster = TRUE,
                        show_row_names = FALSE, show_column_names = FALSE,
                        #column_title_rot = 45,
                        column_title_gp = gpar(fontsize = 0),
                        row_title_rot = 0,
                        row_gap = unit(0.5, "mm"),
                        column_gap = unit(0.5, "mm"),
                        heatmap_legend_param = list(title = '',
                                                    title_gp = gpar(fontface = 1),
                                                    labels_gp = gpar(fontsize = 6)))

pdf("mk_merge.pdf", 
    width = 4.7, height = 3.2)
draw(ht)
dev.off()

ht2 <- signature_heatmap(mk_NN,scale = "sphere", text = FALSE , top = 5,
                         show_row_names = FALSE, show_column_names = FALSE,
                         #column_title_rot = 45,
                         column_title_gp = gpar(fontsize = 0),
                         row_title_rot = 0,
                         row_gap = unit(0.5, "mm"),
                         column_gap = unit(0.5, "mm"),
                         heatmap_legend_param = list(title = '',
                                                     title_gp = gpar(fontface = 1),
                                                     labels_gp = gpar(fontsize = 6)))

pdf("mk_NN.pdf", 
    width = 4.7, height = 3.2)
draw(ht2, padding = unit(c(2, 2, 2, 60), "mm")) 
dev.off()

ht3 <- signature_heatmap(mk_neuron,scale = "sphere", text = FALSE , top = 5,
                         show_row_names = FALSE, show_column_names = FALSE,
                         #column_title_rot = 45,
                         column_title_gp = gpar(fontsize = 0),
                         row_title_rot = 0,
                         row_gap = unit(0.5, "mm"),
                         column_gap = unit(0.5, "mm"),
                         heatmap_legend_param = list(title = '',
                                                     title_gp = gpar(fontface = 1),
                                                     labels_gp = gpar(fontsize = 6)))


pdf("mk_neuron.pdf", 
    width = 4.7, height = 3.2)
draw(ht3)
dev.off()