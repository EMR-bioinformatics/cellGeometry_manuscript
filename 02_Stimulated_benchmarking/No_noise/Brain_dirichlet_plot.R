####Background####

# benchmark plots of brain simulated data

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain")

library(cellGeometry)

cellgeo <- readRDS("cellgeo_dirichlet_output.rds")

music2 <- readRDS("music_dirichlet_output.rds")
music2_NN <- readRDS("music_dirichlet_output_NN.rds")

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
  
  metrics_percent <- metric_set(sim_percent, propdf)
  
  music2[["Times 1"]][[x]]$metrics_percent <- metrics_percent
  
}

music2_metric <- rbind(arrange_metric(music2, "MuSiC"),
                       arrange_metric(music2_NN, "MuSiC"))

metric_percent <- rbind(cellgeo_metric,
                        music2_metric)

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
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 2, scale = "free_y")+
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
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 2, scale = "free_y")+
  labs(x = "", y = "Rsq") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))

ggarrange(Pv2,P2v2, ncol = 1, nrow = 2,
          align = "hv")

# shared y axis

# RMSE
ggplot(data = metric_percent_NN) +
  geom_jitter(aes(x = Celltype, y = RMSEv2, group = Celltype, color = Method),
              alpha = 0.8, 
              width = 0.2)+
  geom_boxplot(aes(x = Celltype, y = RMSE), 
               alpha = 0.5) + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  labs(x = "", y = "RMSE") +
  facet_wrap(~Method, ncol = 2)+
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))


# Rsq
ggplot(data = metric_percent_NN) +
  geom_jitter(aes(x = Celltype, y = Rsq, group = Celltype, color = Method),
              alpha = 0.8, 
              width = 0.2)+
  geom_boxplot(aes(x = Celltype, y = Rsq), 
               alpha = 0.5) + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  labs(x = "", y = "Rsq") +
  facet_wrap(~Method, ncol = 2)+
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))


####neuron metric####

#only three reps for neuron
metric_percent_neuron <- metric_percent[metric_percent$Rep %in% c(1:3), ]
metric_percent_neuron <- metric_percent_neuron[grepl("Human", metric_percent_neuron$Cluster), ]

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
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 2, scale = "free_y")+
  labs(x = "", y = "RMSE") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

P4v2 <- ggplot(metric_percent_neuron, aes(x = Celltype, y = Rsq, 
                                          group = Celltype)) +
  geom_jitter(aes(color = Method), alpha = 0.6, width = 0.1)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 2, scale = "free_y")+
  labs(x = "", y = "Rsq") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

ggarrange(P3v2,P4v2, ncol = 1, nrow = 2,
          align = "hv")


#with removal of outliers/shared y axis
filter_lims <- function(x){
  if(length(x) > 3){
    l <- boxplot.stats(x)$stats[1]
    u <- boxplot.stats(x)$stats[5]
    
    for (i in 1:length(x)){
      x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
    }
    x
  }
  else(
    x
  )
}

#https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

library(dplyr)

metric_percent_neuron <- metric_percent_neuron %>% group_by(Method, Celltype) %>%
  mutate(RMSEv2 = filter_lims(RMSE),
         Rsqv2 = filter_lims(Rsq))

metric_percent_neuron$RMSEv3 <- NA

filt <- max(metric_percent_neuron$RMSEv2, na.rm = TRUE)

metric_percent_neuron$RMSEv3[metric_percent_neuron$RMSE < filt] <- metric_percent_neuron$RMSE[metric_percent_neuron$RMSE < filt]

# RMSE

ggplot(data = metric_percent_neuron) +
  geom_jitter(aes(x = Celltype, y = RMSEv3, group = Celltype, color = Method),
              alpha = 0.8, 
              width = 0.2)+
  stat_summary(aes(x = Celltype, y = RMSE), 
               alpha = 0.5,
               fun.data = calc_boxplot_stat, geom="boxplot") + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  labs(x = "", y = "RMSE") +
  facet_wrap(~Method, ncol = 2)+
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        plot.margin = unit(c(0, 0, 0, 1.5), "cm"))

# Rsq

metric_percent_neuron$Rsqv3 <- NA

filtv2 <- min(metric_percent_neuron$Rsqv2, na.rm = TRUE)

metric_percent_neuron$Rsqv3[metric_percent_neuron$Rsq > filtv2] <- metric_percent_neuron$Rsq[metric_percent_neuron$Rsq > filtv2]

ggplot(data = metric_percent_neuron) +
  geom_jitter(aes(x = Celltype, y = Rsqv3, group = Celltype, color = Method),
              alpha = 0.8, 
              width = 0.2)+
  stat_summary(aes(x = Celltype, y = Rsq), 
               alpha = 0.5,
               fun.data = calc_boxplot_stat, geom="boxplot") + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  labs(x = "", y = "Rsq") +
  facet_wrap(~Method, ncol = 2)+
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        plot.margin = unit(c(0, 0, 0, 1.5), "cm"))

####scatter plots####

plot_setv2 <- function(obs, pred, mfrow = NULL,
                       show_zero = FALSE,
                       show_identity = FALSE,
                       cols = NULL,
                       colour = "blue",
                       title = "", cex.title = 1, ...) {
  if (!identical(dim(obs), dim(pred))) stop("incompatible dimensions")
  if (is.null(cols)) cols <- TRUE
  subclasses <- colnames(obs)[cols]
  nr1 <- ceiling(sqrt(length(subclasses)))
  nr2 <- ceiling(length(subclasses) / nr1)
  if (is.null(mfrow)) mfrow <- c(nr1, nr2)
  oma <- par("oma")
  if (title != "" & oma[3] < 1.5) oma[3] <- 1.5
  op <- par(bty = "l", mgp = c(2.2, 0.6, 0), tcl = -0.3, oma = oma,
            mar = c(3.7, 3.7, 1.5, 1.1), mfrow = mfrow)
  on.exit(par(op))
  scheme <- rev(hue_pal(h = c(0, 270), c = 120)(11))
  xlim <- ylim <- NULL
  new.args <- list(...)
  
  for (i in subclasses) {
    if (is.na(i)) {plot.new(); next}
    if (show_zero) {
      xr <- range(obs[, i], na.rm = TRUE)
      xlim <- c(min(xr[1], 0), xr[2])
      yr <- range(pred[, i], na.rm = TRUE)
      ylim <- c(min(yr[1], 0), yr[2])
    }
    args <- list(x = obs[, i], y = pred[, i], cex = 0.8, pch = 16, las = 1,
                 xlab = i %>% str_wrap(width = 20), ylab = "Predicted", xlim = xlim, ylim = ylim)
    if (length(new.args)) args[names(new.args)] <- new.args
    do.call(plot, args)
    # fit <- lm(pred[, i] ~ obs[, i])
    # rsq <- summary(fit)$r.squared
    Rsq <- Rsq(obs[ , i], pred[, i])
    col <- if (colour == "rainbow") scheme[ceiling(rsq*10) +1] else colour
    #abline(fit, col = col, lwd = 1.5)
    if (show_identity) abline(0, 1, col = "grey50", lty = 2)
    mtext(bquote(R^2 == .(format(Rsq, digits = 3))), cex = 0.8, adj = 0.04)
  }
  mtext(title, outer = TRUE, cex = cex.title * par("cex"), adj = 0.05, line = 0)
}

Rsq <- function(obs, pred) {
  rss <- sum((pred - obs)^2)
  tss <- sum((obs - mean(obs))^2)
  1 - rss/tss
}

#####neuronal#####

sim_percent <- sim_percent_merge[ , grepl("Human", colnames(sim_percent_merge))]
cellgeo_percent <- as.data.frame(cellgeo$`Times 1`$`Rep 1`$output$subclass$percent)

plot_setv2(sim_percent[ , 1:35], cellgeo_percent[ , 1:35],
           show_identity = TRUE, show_zero = TRUE, mfrow = c(7,5))


plot_setv2(sim_percent[ , 36:70], cellgeo_percent[ , 36:70],
           show_identity = TRUE, show_zero = TRUE,mfrow = c(7,5))


plot_setv2(sim_percent[ , 71:105], cellgeo_percent[ , 71:105],
           show_identity = TRUE, show_zero = TRUE,mfrow = c(7,5))


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

plot_setv2(sim_percent[ , 1:35], music_percent[ , 1:35],
           show_identity = TRUE, show_zero = TRUE,mfrow = c(7,5))


plot_setv2(sim_percent[ , 36:70], music_percent[ , 36:70],
           show_identity = TRUE, show_zero = TRUE,mfrow = c(7,5))

plot_setv2(sim_percent[ , 71:105], music_percent[ , 71:105],
           show_identity = TRUE,show_zero = TRUE, mfrow = c(7,5))


#####non-neuronal#####

sim_percentNN <- sim_percent_merge[ , !grepl("Human", colnames(sim_percent_merge))]

all(colnames(sim_percentNN) == colnames(cellgeo_percent)[106:116])
#TRUE

plot_setv2(sim_percentNN, cellgeo_percent[ , 106:116],
           show_identity = TRUE, show_zero = TRUE,mfrow = c(2,6))


prop_arrangev2 <- function(df){
  common <- intersect(colnames(sim_percentNN), colnames(df))
  propdf <- df[ , colnames(df) %in% common]
  
  missing <- setdiff(colnames(sim_percentNN), colnames(propdf))
  
  propdf <- as.data.frame(propdf)
  
  if(length(missing) > 0){
    for(i in missing){
      propdf[ , i] <- 0
    }
  }
  
  propdf <- propdf[ , match(colnames(sim_percentNN), colnames(propdf))]
  propdf
}

music_percentNN <- prop_arrangev2(music2_NN$`Times 1`$`Rep 1`$output$Est.prop.weighted * 100)

all(colnames(sim_percentNN) == colnames(music_percentNN))
#TRUE

plot_setv2(sim_percentNN, music_percentNN,
           show_identity = TRUE, show_zero = TRUE,mfrow = c(2,6))

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




