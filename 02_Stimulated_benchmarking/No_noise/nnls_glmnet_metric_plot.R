# attempt to combine nnls and glmnet metrics with NMF 

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/NNLS")

####organise data####

path <- "/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/Brain/metric_plots/NMF"

NMF_metrics <- readRDS(paste0(path, "/NMF_metrics.rds"))

arrange_metric <- function(list, method, data){
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
  df$Data <- data
  df
}

amp_nnls <- readRDS("nnls_AMP.rds")
typist_nnls <- readRDS("nnls_typist.rds")
tabula_nnls <- readRDS("nnls_tabula.rds")
brain_nnls <- readRDS("nnls_brain.rds")

library(tibble)

nnls_metrics <- rbind(arrange_metric(amp_nnls, "NNLS", "AMP"),
                      arrange_metric(typist_nnls, "NNLS", "Cell typist"),
                      arrange_metric(tabula_nnls, "NNLS", "Tabula Sapiens"),
                      arrange_metric(brain_nnls, "NNLS", "Human Brain\nCell Atlas"))

nnls_metrics$Data_info <- NA

#glmnet using subsetted cellGeometry signature to be comparable like nnls
amp_glmnet <- readRDS("glmnet_amp.rds")
typist_glmnet <- readRDS("glmnet_typist.rds")
tabula_glmnet <- readRDS("glmnet_tabula.rds")
brain_glmnet <- readRDS("glmnet_brain.rds")

arrange_metricv2 <- function(list, method, data){
  df <- data.frame()
  for(i in names(list)){
    for(j in names(list[[i]])){
      for(x in names(list[[i]][[j]])){
      temp <- as.data.frame(list[[i]][[j]][[x]]$metrics_percent)
      temp <- add_column(temp,
                         "Cluster" = rownames(temp),
                         .before = 1)
      temp$Times <- as.numeric(gsub("Times ", "", i))
      temp$Rep <- as.numeric(gsub("Rep ", "", x))
      temp$Method <- paste0(method, " \u03b1", "=", gsub("alpha = ", "", j))
      df <- rbind(df,
                  temp)
      }
    }
  }
  
  df$Data <- data
  df
}

glmnet_metrics <- rbind(arrange_metricv2(amp_glmnet, "glmnet", "AMP"),
                      arrange_metricv2(typist_glmnet, "glmnet", "Cell typist"),
                      arrange_metricv2(tabula_glmnet, "glmnet", "Tabula Sapiens"),
                      arrange_metricv2(brain_glmnet, "glmnet", "Human Brain\nCell Atlas"))

glmnet_metrics$Data_info <- NA

metric_all <- rbind(NMF_metrics,
                    nnls_metrics,
                    glmnet_metrics[glmnet_metrics$Method %in% c("glmnet α=0",
                                                                "glmnet α=1"), ])

# saveRDS(metric_all, "comb_metrics.rds")

####plot####

library(ggh4x)

ggplot(data = metric_all) +
  geom_jitter(aes(x = Method, y = RMSE, color = Method), alpha = 0.2, 
              width = 0.2)+
  geom_boxplot(aes(x = Method, y = RMSE), 
               alpha = 0.5, outlier.shape = NA) + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "Brunet" = "grey60",
                                "KL" = "burlywood3",
                                "Lee" = "slategray3",
                                "RcppML" = "thistle3",
                                "NNLS" = "darkseagreen3",
                                "glmnet α=0" = "olivedrab2",
                                "glmnet α=1" = "olivedrab4"),
                     guide = "none") +
  labs(x = "", y = "RMSE") +
  coord_cartesian(clip = "off") +
  facet_grid2(~Data, scales = "free", space = "free_x", independent = "y")+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #strip.text = element_text(size = 12),
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.3, "cm"),
        panel.spacing.x = unit(1, "cm"))
#800x380

ggplot(data = metric_all) +
  geom_jitter(aes(x = Method, y = Rsq, color = Method), alpha = 0.2, 
              width = 0.2)+
  geom_boxplot(aes(x = Method, y = Rsq), 
               alpha = 0.5, outlier.shape = NA) + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "Brunet" = "grey60",
                                "KL" = "burlywood3",
                                "Lee" = "slategray3",
                                "RcppML" = "thistle3",
                                "NNLS" = "darkseagreen3",
                                "glmnet α=0" = "olivedrab2",
                                "glmnet α=1" = "olivedrab4"),
                     guide = "none") +
  labs(x = "", y = "Rsq") +
  coord_cartesian(clip = "off") +
  facet_grid2(~Data, scales = "free", space = "free_x", independent = "y")+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #strip.text = element_text(size = 12),
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.3, "cm"),
        panel.spacing.x = unit(1, "cm"))
#800x380

####cropped plots####

metric_plot <- function(data, metric, ylim){
  ggplot(data = metric_all[metric_all$Data == data, ]) +
    geom_jitter(aes(x = Method, y = .data[[metric]], color = Method), alpha = 0.2, 
                width = 0.2)+
    geom_boxplot(aes(x = Method, y = .data[[metric]]), 
                 alpha = 0.5, outlier.shape = NA) + 
    scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                  "Brunet" = "grey60",
                                  "KL" = "burlywood3",
                                  "Lee" = "slategray3",
                                  "RcppML" = "thistle3",
                                  "NNLS" = "darkseagreen3",
                                  "glmnet α=0" = "olivedrab2",
                                  "glmnet α=1" = "olivedrab4"),
                       guide = "none") +
    labs(x = "", y = metric) +
    coord_cartesian(ylim = ylim) +
    theme_classic()+
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(0.1, 0.1, 0.1, 0.5, "cm"),
          panel.spacing.x = unit(1, "cm"))
}

library(ggpubr)

ggarrange(metric_plot("AMP", "RMSE", range(metric_all$RMSE[metric_all$Data == "AMP"])),
          metric_plot("Cell typist", "RMSE", range(metric_all$RMSE[metric_all$Data == "Cell typist"])) + ylab(""),
          metric_plot("Tabula Sapiens", "RMSE", range(0, 2)) + ylab(""),
          metric_plot("Human Brain\nCell Atlas", "RMSE", range(0, 4)) + ylab(""),
          ncol = 4)
#800x380

min(metric_all$Rsq[metric_all$Data == "Tabula Sapiens" & metric_all$RMSE < 2])
#-45.2424

min(metric_all$Rsq[metric_all$Data == "Human Brain\nCell Atlas" & metric_all$RMSE < 2])
# -18.57613

ggarrange(metric_plot("AMP", "Rsq", range(metric_all$Rsq[metric_all$Data == "AMP"])),
          metric_plot("Cell typist", "Rsq", range(metric_all$Rsq[metric_all$Data == "Cell typist"])) + ylab(""),
          metric_plot("Tabula Sapiens", "Rsq", c(-46, 1)) + ylab(""),
          metric_plot("Human Brain\nCell Atlas", "Rsq", c(-19, 1)) + ylab(""),
          ncol = 4)
#800x380

####alpha####

get_rsq <- function(list){
  df <- data.frame()
  for(i in names(list)){
    for(j in names(list[[i]])){
      for(x in names(list[[i]][[j]])){
        temp <- data.frame("Rsq_percent" = list[[i]][[j]][[x]]$rsq_percent,
                           "Rsq_counts" = list[[i]][[j]][[x]]$rsq_counts,
                           "Times" = as.numeric(gsub("Times ", "", i)),
                           "Rep" = as.numeric(gsub("Rep ", "", x)),
                           "Alpha" = gsub("alpha = ", "", j))
        df <- rbind(df,
                    temp)
      }
    }
  }
  
  df
}

typist_rsq <- get_rsq(typist_glmnet)

library(dplyr)

typist_rsq_sem <- typist_rsq %>%
  group_by(Alpha) %>%
  summarise("mean_percent" = mean(Rsq_percent),
            "mean_counts" = mean(Rsq_counts),
            "SD_percent" = sd(Rsq_percent),
            "SD_counts" = sd(Rsq_counts),
            "SEM_percent" = SD_percent/sqrt(n()),
            "SEM_counts" = SD_counts/sqrt(n()))

ggplot(typist_rsq_sem, aes(x = Alpha, y = mean_percent)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_percent - SEM_percent, ymax = mean_percent + SEM_percent),
                width = 0.2) +
  labs(x = "\u03b1", y = "Rsq") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
#400x250


