####Background####

#benchmark plots of Tabula Sapiens simulated data (now with DWLS)

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2")

library(cellGeometry)
cellgeo <- readRDS("cellgeo_dirichlet_output_500.rds")
music2 <- readRDS("music_dirichlet_output.rds")
DWLS <- readRDS("DWLS_dirichlet_output.rds")

mk <- readRDS("tabula_markers_dualmeans_rerun.rds")

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Tabula/workstation2/tabula_simulated_workstation2.rdata")

celltypedf <- readRDS("DWLS_celltype_convert.rds")

#get metrics for DWLS 

metric_out <- function(list){
  for(i in names(list)){
    for(x in names(list[[i]])){
      df <- t(list[[i]][[x]]$output * 100)
      colnames(df) <- celltypedf$Orig[match(colnames(df),
                                            celltypedf$convert)]
      df <- df[ , colnames(sim_percent)]
      
      list[[i]][[x]]$metrics_percent <- metric_set(sim_percent, df)
    }
  }
  list
}

DWLS <- metric_out(DWLS)

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

DWLS_metric <- arrange_metric(DWLS, "DWLS")

metric_percent <- rbind(cellgeo_metric,
                        music2_metric,
                        DWLS_metric)

metric_percent$Method <- factor(metric_percent$Method,
                                levels = c("CellGeometry",
                                           "MuSiC",
                                           "DWLS"))

metric_percent$Data <- "Tabula"

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
                                "MuSiC" = "#00BF7D",
                                "DWLS" = "#e86bf3"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 3, scale = "free_y")+
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
                                "MuSiC" = "#00BF7D",
                                "DWLS" = "#e86bf3"),
                     guide = "none") +
  facet_wrap(~Method, ncol = 3, scale = "free_y")+
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

ggplot(data = metric_percent[metric_percent$Method != "DWLS", ]) +
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
                                "MuSiC" = "#00BF7D",
                                "DWLS" = "#e86bf3"),
                     guide = "none") +
  labs(x = "", y = "Rsq") +
  facet_wrap2(~Method, ncol = 3, axes = "all")+
  theme_classic()+
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))
#panel.spacing.x = unit(1.5, "cm"))

####scatter plot####

plot_pred(sim_percent, cellgeo$`Times 3`$`Rep 1`$output$subclass$percent, mk) + labs(x = "Observed (%)", y = "Predicted (%)")

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

plot_pred(sim_percent, as.matrix(music_percent), mk, ellipse = 3) + labs(x = "Observed (%)", y = "Predicted (%)")

prop_arrangev2 <- function(df){
  df <- t(df * 100)
  colnames(df) <- celltypedf$Orig[match(colnames(df),
                                        celltypedf$convert)]
  df <- df[ , colnames(sim_percent)]
}

DWLS_percent <- prop_arrangev2(DWLS$`Times 3`$`Rep 1`$output)

all(colnames(sim_percent) == colnames(DWLS_percent))
#TRUE

plot_pred(sim_percent, as.matrix(DWLS_percent), mk) + labs(x = "Observed (%)", y = "Predicted (%)")
