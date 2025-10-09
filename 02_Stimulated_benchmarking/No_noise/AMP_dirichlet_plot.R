####Background####

# benchmark plots of AMP simulated data

library(cellGeometry)

cellgeo <- readRDS("cellgeo_dirichlet_output500.rds")

music2 <- readRDS("Music_dirichlet_output.rds")
DWLS <- readRDS("DWLS_dirichlet_output.rds")

Lin <- readRDS("Lin_dirichlet_output.rds")

col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))

cellgeo_percent <- as.data.frame(cellgeo$`Times 30`$`Rep 2`$output$subclass$percent)
music_percent <- as.data.frame(music2$`Times 30`$`Rep 2`$output$Est.prop.weighted * 100)
DWLS_percent <- as.data.frame(t(DWLS$`Times 30`$`Rep 2`$output  * 100))
Lin_percent <- as.data.frame(Lin$`Times 30`$`Rep 2`$output * 100)
colnames(Lin_percent) <- gsub("\\.", "-", colnames(Lin_percent))
rownames(Lin_percent) <- gsub("\\.", "-", rownames(Lin_percent))
load("simulated_dirichlet.rdata")

library(tidyverse)

arrange_data <- function(df, method){
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = "Percentage"))
  test$Method <- method
  test
}

comb <- rbind(arrange_data(as.data.frame(sim_percent), "Ground truth"),
              arrange_data(cellgeo_percent, "CellGeometry"),
              arrange_data(music_percent, "MuSiC"),
              arrange_data(DWLS_percent, "DWLS"),
              arrange_data(Lin_percent, "LinDeconSeq"))

comb$Cluster <- factor(comb$Cluster,
                       levels = col_order)

comb$Method <- factor(comb$Method,
                      levels = c("Ground truth",
                                 "CellGeometry",
                                 "MuSiC",
                                 "DWLS",
                                 "LinDeconSeq"))

comb$Celltype <- NA
comb$Celltype[grepl("F", comb$Cluster)] <- "Fibroblast"
comb$Celltype[grepl("M", comb$Cluster)] <- "Macrophage"
comb$Celltype[grepl("B", comb$Cluster)] <- "B cell"
comb$Celltype[grepl("T", comb$Cluster)] <- "T cell"

comb$Celltype <- factor(comb$Celltype,
                        levels = c("Fibroblast", "Macrophage",
                                   "T cell", "B cell"))

####ridge plot####

library(ggridges)
library(scales)
library(ggplot2)

amp_scheme <- c("SC-F1" = "palegreen1",
                "SC-F2" = "seagreen3",
                "SC-F3" = "green4",
                "SC-F4" = "darkgreen",
                "SC-M1" = "rosybrown1",
                "SC-M2" = "tomato1",
                "SC-M3" = "firebrick1",
                "SC-M4" = "firebrick3",
                "SC-T1" = "plum1",
                "SC-T2" = "plum3",
                "SC-T3" = "mediumorchid3",
                "SC-T4" = "purple2",
                "SC-T5" = "darkorchid1",
                "SC-T6" = "darkmagenta",
                "SC-B1" = "lightblue",
                "SC-B2" = "skyblue1",
                "SC-B3" = "royalblue1",
                "SC-B4" = "darkblue")

ggplot(comb,
       aes(x = Percentage, y = Cluster,
           fill = Cluster)) +
  geom_density_ridges(alpha = 0.6, na.rm = TRUE) +
  scale_fill_manual(values = amp_scheme)+
  theme_ridges() + 
  scale_x_continuous(breaks = c(0, 10, 20, 30))+
  labs(x = "Cell (%)", y = "")+
  facet_wrap(~Method, ncol = 5, scales = "free_x")+
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        strip.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text = element_text(size = 12, angle = 45, 
                                  hjust = 0.5, vjust = 0.5))

####metric boxplots####

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
  df$Cluster <- factor(df$Cluster,
                       levels = col_order)
  df
}

library(tibble)

cellgeo_metric <- arrange_metric(cellgeo, "CellGeometry")

music2_metric <- arrange_metric(music2, "MuSiC")

DWLS_metric <- arrange_metric(DWLS, "DWLS")

Lin_metric <-arrange_metric(Lin, "LinDeconSeq")

metric_percent <- rbind(cellgeo_metric,
                        music2_metric,
                        DWLS_metric,
                        Lin_metric)

metric_percent$Celltype <- NA
metric_percent$Celltype[grepl("F", metric_percent$Cluster)] <- "Fibroblast"
metric_percent$Celltype[grepl("M", metric_percent$Cluster)] <- "Macrophage"
metric_percent$Celltype[grepl("B", metric_percent$Cluster)] <- "B cell"
metric_percent$Celltype[grepl("T", metric_percent$Cluster)] <- "T cell"

metric_percent$Celltype <- factor(metric_percent$Celltype,
                                  levels = c("Fibroblast", "Macrophage",
                                             "T cell", "B cell"))

metric_percent$Method <- factor(metric_percent$Method,
                                c("CellGeometry", "MuSiC", "DWLS", "LinDeconSeq"))

library(ggh4x)

# RMSE boxplot - cropped y axis

ggplot(metric_percent, aes(x = Celltype, y = RMSE,
                           color = Cluster, group = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7, width = 0.1)+
  scale_color_manual(values = amp_scheme, name = "", guide = "none")+
  facet_wrap2(~Method, ncol = 4, axes = "all")+
  coord_cartesian(ylim = c(0, 6.5))+
  labs(x = "", y = "RMSE") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.x = element_blank(),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        panel.spacing.x = unit(1.5, "cm"))

# RMSE boxplot - independent full y axis

Pv2 <- ggplot(metric_percent, aes(x = Celltype, y = RMSE,
                           color = Cluster, group = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7, width = 0.1)+
  scale_color_manual(values = amp_scheme, name = "", guide = "none")+
  facet_wrap(~Method, ncol = 4, scale = "free_y")+
  labs(x = "", y = "RMSE") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))

# Rsq boxplot - cropped y axis 

ggplot(metric_percent, aes(x = Celltype, y = Rsq,
                           color = Cluster, group = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7, width = 0.1)+
  scale_color_manual(values = amp_scheme, name = "", guide = "none")+
  facet_wrap2(~Method, ncol = 4, axes = "all")+
  coord_cartesian(ylim = c(-4 ,1))+
  labs(x = "", y = "Rsq") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        panel.spacing.x = unit(1.5, "cm"))

# Rsq boxplot - independent full y axis 

P2v2 <- ggplot(metric_percent, aes(x = Celltype, y = Rsq,
                           color = Cluster, group = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7, width = 0.1)+
  scale_color_manual(values = amp_scheme, name = "", guide = "none")+
  facet_wrap(~Method, ncol = 4, scales = "free_y")+
  labs(x = "", y = "Rsq") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))

library(ggpubr)

ggarrange(Pv2,P2v2, ncol = 1, nrow = 2,
          align = "hv")

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
                 xlab = i, ylab = "Predicted", xlim = xlim, ylim = ylim)
    if (length(new.args)) args[names(new.args)] <- new.args
    do.call(plot, args)
    Rsq <- Rsq(obs[ , i], pred[, i])
    col <- if (colour == "rainbow") scheme[ceiling(rsq*10) +1] else colour
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

plot_setv2(sim_percent, cellgeo_percent, show_identity = TRUE, 
           show_zero = TRUE, mfrow = c(3, 6))

plot_setv2(sim_percent, music_percent, show_identity = TRUE, 
           show_zero = TRUE, mfrow = c(3, 6))

plot_setv2(sim_percent, DWLS_percent, show_identity = TRUE, 
           show_zero = TRUE,mfrow = c(3, 6))
#750x350

plot_setv2(sim_percent, Lin_percent, show_identity = TRUE,
           show_zero = TRUE, mfrow = c(3, 6))
#750x350






