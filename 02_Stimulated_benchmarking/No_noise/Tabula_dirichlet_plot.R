####Background####

#benchmark plots of Tabula Sapiens simulated data

cellgeo <- readRDS("cellgeo_dirichlet_output_500.rds")
music2 <- readRDS("music_dirichlet_output.rds")
mk <- readRDS("tabula_markers_dualmeans_rerun.rds")

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

library(tibble)

cellgeo_metric <- arrange_metric(cellgeo, "CellGeometry")

music2_metric <- arrange_metric(music2, "MuSiC")

metric_percent <- rbind(cellgeo_metric,
                        music2_metric)

####metric boxplots####

cell_groups <- as.data.frame(mk$cell_table)
colnames(cell_groups) <- "Group"

metric_percent$Celltype <- cell_groups$Group[match(metric_percent$Cluster,
                                                   rownames(cell_groups))]

library(ggplot2)

# full RMSE boxplot 

Pv2 <- ggplot(metric_percent, aes(x = Celltype, y = RMSE, group = Celltype)) +
  geom_jitter(aes(color = Method), alpha = 0.15, width = 0.1, size = 0.8)+
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

# full Rsq boxplot 

P2v2 <- ggplot(metric_percent, aes(x = Celltype, y = Rsq, group = Celltype)) +
  geom_jitter(aes(color = Method), alpha = 0.15, width = 0.1, size = 0.8)+
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

library(ggpubr)

ggarrange(Pv2,P2v2, ncol = 1, nrow = 2,
          align = "hv")

# cropped y axis (outliers removed)

filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

library(dplyr)

metric_percent <- metric_percent %>% group_by(Method, Celltype) %>%
  mutate(RMSEv2 = filter_lims(RMSE),
         Rsqv2 = filter_lims(Rsq))

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

metric_percent$RMSEv3 <- NA

filt <-  max(metric_percent$RMSEv2, na.rm = TRUE)

metric_percent$RMSEv3[metric_percent$RMSE < filt] <- metric_percent$RMSE[metric_percent$RMSE < filt]

metric_percent$Rsqv3 <- NA

filtv2 <-  min(metric_percent$Rsqv2, na.rm = TRUE)

metric_percent$Rsqv3[metric_percent$Rsq > filtv2] <- metric_percent$Rsq[metric_percent$Rsq > filtv2]

library(ggh4x)

# RMSE

ggplot(data = metric_percent) +
  geom_jitter(aes(x = Celltype, y = RMSEv3, group = Celltype, color = Method),
              alpha = 0.2, 
              width = 0.1)+
  stat_summary(aes(x = Celltype, y = RMSE), 
               alpha = 0.5,
               fun.data = calc_boxplot_stat, geom="boxplot") + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  facet_wrap2(~Method, ncol = 2, axes = "all")+
  labs(x = "", y = "RMSE") +
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        panel.spacing.x = unit(1.5, "cm"))

# Rsq

ggplot(data = metric_percent) +
  geom_jitter(aes(x = Celltype, y = Rsqv3, group = Celltype, color = Method),
              alpha = 0.2, 
              width = 0.1)+
  stat_summary(aes(x = Celltype, y = Rsq), 
               alpha = 0.5,
               fun.data = calc_boxplot_stat, geom="boxplot") + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D"),
                     guide = "none") +
  labs(x = "", y = "Rsq") +
  facet_wrap2(~Method, ncol = 2, axes = "all")+
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        panel.spacing.x = unit(1.5, "cm"))

####scatter plot####

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
                 xlab = i %>% str_wrap(width = 25), # was 30
                 ylab = "Predicted", xlim = xlim, ylim = ylim)
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

load("tabula_simulated_workstation2.rdata")

cellgeo_percent <- as.data.frame(cellgeo$`Times 3`$`Rep 1`$output$subclass$percent)

all(colnames(sim_percent) == colnames(cellgeo_percent))
#TRUE

library(stringr)

plot_setv2(sim_percent[ , 1:45], cellgeo_percent[ , 1:45], show_zero = TRUE,
           show_identity = TRUE, mfrow = c(9,5))


plot_setv2(sim_percent[ , 46:90], cellgeo_percent[ , 46:90], show_zero = TRUE,
           show_identity = TRUE, mfrow = c(9,5))


plot_setv2(sim_percent[ , 91:135], cellgeo_percent[ , 91:135], show_zero = TRUE,
           show_identity = TRUE, mfrow = c(9,5))


plot_setv2(sim_percent[ , 136:171], cellgeo_percent[ , 136:171], show_zero = TRUE,
           show_identity = TRUE, mfrow = c(9,5))


prop_arrange <- function(df){
  common <- intersect(colnames(sim_percent), colnames(df))
  propdf <- df[ , colnames(df) %in% common]
  
  missing <- setdiff(colnames(sim_percent), colnames(propdf))
  print(missing)
  propdf <- as.data.frame(propdf)
  
  if(length(missing) > 0){
    for(i in missing){
      propdf[ , i] <- 0
    }
  }
  
  propdf <- propdf[ , match(colnames(sim_percent), colnames(propdf))]
  propdf
}

music_percent <- prop_arrange(music2$`Times 3`$`Rep 1`$output$Est.prop.weighted * 100)

all(colnames(sim_percent) == colnames(music_percent))
#TRUE

plot_setv2(sim_percent[ , 1:45], music_percent[ , 1:45], show_zero = TRUE,
           show_identity = TRUE, mfrow = c(9,5))


plot_setv2(sim_percent[ , 46:90], music_percent[ , 46:90], show_zero = TRUE,
           show_identity = TRUE, mfrow = c(9,5))


plot_setv2(sim_percent[ , 91:135], music_percent[ , 91:135], show_zero = TRUE,
           show_identity = TRUE, mfrow = c(9,5))


plot_setv2(sim_percent[ , 136:171], music_percent[ , 136:171], show_zero = TRUE,
           show_identity = TRUE, mfrow = c(9,5))



