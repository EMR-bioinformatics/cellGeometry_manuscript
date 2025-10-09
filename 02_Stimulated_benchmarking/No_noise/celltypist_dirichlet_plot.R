####Background####

#benchmark plots of Cell Typist blood simulated data

library(cellGeometry)

cellgeo <- readRDS("cellgeo_dirichlet_output.rds")
music2 <- readRDS("music2_output.rds")
DWLS <- readRDS("DWLS_dirichlet_output.rds")
Lin <- readRDS("LinDeconSeq_output.rds")

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
music2_metric <- arrange_metric(music2, "MuSiC")
DWLS_metric <- arrange_metric(DWLS, "DWLS")
Lin_metric <-arrange_metric(Lin, "LinDeconSeq")

metric_percent <- rbind(cellgeo_metric,
                        music2_metric,
                        DWLS_metric,
                        Lin_metric)

col_order <- c(
  "Memory B cells", "Naive B cells", 
  "Plasma cells",
  "DC", "DC2",
  "pDC", "GMP",
  "Mast cells",
  "Classical monocytes", "Non-classical monocytes",
  "Macrophages",
  "Follicular helper T cells", "MAIT cells", "Regulatory T cells",
  "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells",
  "Tem/Effector cytotoxic T cells", "Tem/Effector helper T cells",
  "Cycling T cells",
  "CD16+ NK cells", "CD16- NK cells", "NK cells", "ILC",
  "Early MK",
  "HSC/MPP")

metric_percent$Cluster <- factor(metric_percent$Cluster,
                                 levels = col_order)

load("simulated_dirichlet.rdata")

####Ridge plot####

cellgeo_percent <- as.data.frame(cellgeo$`Times 30`$`Rep 1`$output$subclass$percent)

prop_arrange <- function(df){
  prop <- df$Est.prop.weighted * 100
  common <- intersect(colnames(sim_percent), colnames(prop))
  
  prop <- prop[ , colnames(prop) %in% common]
  
  missing <- setdiff(colnames(sim_percent), colnames(prop))
  
  propdf <- as.data.frame(prop)
  
  if(length(missing) > 0){
    for(i in missing){
      propdf[ , i] <- 0
    }
  }
  
  propdf[ , match(colnames(sim_percent), colnames(propdf))]
}


music_percent <- as.data.frame(prop_arrange(music2$`Times 30`$`Rep 1`$output))

prop_arrangev2 <- function(df){
  colnames(df) <- gsub("_", " ", colnames(df))
  colnames(df) <- gsub("neg", "-", colnames(df))
  colnames(df)[colnames(df) == "Tem Effector helper T cells"] <- "Tem/Effector helper T cells"
  colnames(df)[colnames(df) == "Nonclassical monocytes"] <- "Non-classical monocytes"
  colnames(df)[colnames(df) == "CD16plus NK cells"] <- "CD16+ NK cells" 
  colnames(df)[colnames(df) == "Tcm Naive helper T cells"] <- "Tcm/Naive helper T cells" 
  colnames(df)[colnames(df) == "Tcm Naive cytotoxic T cells"] <- "Tcm/Naive cytotoxic T cells" 
  colnames(df)[colnames(df) == "Tem Effector cytotoxic T cells"] <- "Tem/Effector cytotoxic T cells" 
  colnames(df)[colnames(df) == "HSC MPP"] <- "HSC/MPP" 
  
  common <- intersect(colnames(sim_percent), colnames(df))
  df <- df[ , colnames(df) %in% common]
  df[ , match(colnames(sim_percent), colnames(df))]
}


DWLS_percent <- as.data.frame(prop_arrangev2(t(DWLS$`Times 30`$`Rep 1`$output * 100)))
Lin_percent <- as.data.frame(prop_arrangev2(Lin$`Times 30`$`Rep 1`$output * 100))

col_scheme <- c(
  "Memory B" = "skyblue1",
  "Naive B" = "lightblue",
  "Plasma" = "royalblue4",
  "DC" = "sienna4",
  "DC2" = "sienna3",
  "pDC" = "sienna1",
  "GMP" = "orange",
  "Mast" = "brown4",
  "Classical monocytes" = "firebrick3",
  "Non-classical monocytes" = "firebrick1",
  "Macrophages" = "rosybrown1",
  "Follicular helper T" = "mediumorchid3",
  "MAIT" = "violet",
  "Regulatory T" = "plum3",
  "Tcm/Naive cytotoxic T" = "plum1",
  "Tcm/Naive helper T" = "purple2",
  "Tem/Effector cytotoxic T" = "darkorchid1",
  "Tem/Effector helper T" = "darkmagenta",
  "Cycling T" = "#490092",
  "CD16+ NK" = "aquamarine",
  "CD16- NK" = "aquamarine3",
  "NK" = "aquamarine4",
  "ILC" = "lightseagreen",
  "Early MK" = "gray75",
  "HSC/MPP" = "gray50")

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


comb <- comb[comb$Cluster %in% colnames(sim_percent), ]

comb$Method <- factor(comb$Method,
                      levels = c("Ground truth",
                                 "CellGeometry",
                                 "MuSiC",
                                 "DWLS",
                                 "LinDeconSeq"))

library(ggridges)
library(scales)
library(ggplot2)

comb$Cluster <- factor(comb$Cluster,
                       levels = rev(col_order))

comb$Clusterv2 <- as.vector(comb$Cluster)
comb$Clusterv2 <- gsub(" cells", "", comb$Clusterv2)

comb$Clusterv2 <- factor(comb$Clusterv2,
                         levels = rev(gsub(" cells", "", col_order)))

plot <- ggplot(comb,
               aes(x = Percentage, y = Clusterv2,
                   fill = Clusterv2)) +
  geom_density_ridges(alpha = 0.6, na.rm = TRUE) +
  scale_fill_manual(values = col_scheme)+
  theme_ridges() + 
  scale_x_continuous(breaks = c(0, 10, 20, 30))+
  labs(x = "Cell (%)", y = "")+
  facet_wrap(~Method, ncol = 5, scales = "free_x")+
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        strip.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text = element_text(size = 12, angle = 45, hjust = 0.5, vjust = 0.5))


pg <- ggplotGrob(plot)

for(i in which(grepl("strip-t", pg$layout$name))){
  pg$grobs[[i]]$layout$clip <- "off"
}

grid::grid.draw(pg)

####metric boxplots####

mk <- readRDS("typist_cellmarkers_cycling.rds")
cell_groups <- as.data.frame(mk$cell_table)
colnames(cell_groups) <- "Group"
cell_groups$Groupv2 <- as.vector(cell_groups$Group)
cell_groups$Groupv2[cell_groups$Groupv2 == "Cycling cells"] <- "T cells"
cell_groups$Groupv2[cell_groups$Groupv2 %in% c("Monocytes", "DC", "Macrophages", "Mast cells", "pDC")] <- "Myeloid" 
cell_groups$Groupv2[cell_groups$Groupv2 == "Plasma cells"] <- "B cells"
cell_groups$Groupv2[cell_groups$Groupv2 %in% c("Early MK", "HSC/MPP")] <- "Progenitor/misc"
cell_groups$Groupv2[rownames(cell_groups) == "GMP"] <- "Myeloid"

metric_percent$Celltype <- cell_groups$Groupv2[match(metric_percent$Cluster,
                                                     rownames(cell_groups))]

metric_percent$Celltype <- factor(metric_percent$Celltype,
                                  levels = c("B cells", "Myeloid", "T cells", "ILC", "Progenitor/misc"))

metric_percent$Clusterv2 <- as.vector(metric_percent$Cluster)
metric_percent$Clusterv2  <- gsub(" cells", "", metric_percent$Clusterv2 )

metric_percent$Clusterv2  <- factor(metric_percent$Clusterv2 ,
                                    levels = gsub(" cells", "", col_order))

metric_percent$Method <- factor(metric_percent$Method,
                                levels = c("CellGeometry",
                                           "MuSiC",
                                           "DWLS",
                                           "LinDeconSeq"))

# RMSE - fixed y axis 

ggplot(metric_percent, aes(x = Celltype, y = RMSE,
                           color = Clusterv2, group = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7, width = 0.1)+
  scale_color_manual(values = col_scheme, name = "", guide = "none")+
  facet_wrap(~Method, ncol = 4)+
  labs(x = "", y = "RMSE") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14))

# RMSE - independent y axis

Pv2 <- ggplot(metric_percent, aes(x = Celltype, y = RMSE,
                                  color = Clusterv2, group = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7, width = 0.1)+
  scale_color_manual(values = col_scheme, name = "", guide = "none")+
  facet_wrap(~Method, ncol = 4, scale = "free_y")+
  labs(x = "", y = "RMSE") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))

# Rsq - fixed y axis

ggplot(metric_percent, aes(x = Celltype, y = Rsq,
                           color = Clusterv2, group = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7, width = 0.1)+
  scale_color_manual(values = col_scheme, name = "", guide = "none")+
  facet_wrap(~Method, ncol = 4)+
  labs(x = "", y = "Rsq") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14))

# Rsq - independent y axis

P2v2 <- ggplot(metric_percent, aes(x = Celltype, y = Rsq,
                                   color = Clusterv2, group = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7, width = 0.1)+
  scale_color_manual(values = col_scheme, name = "", guide = "none")+
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
                 xlab = gsub(" cells", "", i), ylab = "Predicted", xlim = xlim, ylim = ylim)
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

plot_setv2(sim_percent, cellgeo_percent, show_identity = TRUE, show_zero = TRUE, mfrow = c(5, 5))

plot_setv2(sim_percent, music_percent, show_identity = TRUE, show_zero = TRUE, mfrow = c(5, 5))

plot_setv2(sim_percent, DWLS_percent, show_identity = TRUE, show_zero = TRUE, mfrow = c(5, 5))

plot_setv2(sim_percent, Lin_percent, show_identity = TRUE,show_zero = TRUE,  mfrow = c(5, 5))
