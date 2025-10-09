####Background####

# summarise metric plots

metric_percent <- readRDS("brain_metric_percent.rds")

typist <- readRDS("../Cell_typist/workstation2/Dirichlet/typist_metric_percent.rds")
typist$Data <- "Cell typist"

AMP <- readRDS("../AMP/Dirichlet/AMP_metric_percent_withexp.rds")
AMP$Data <- "AMP"

tabula <- readRDS("../Tabula/workstation2/tabula_metric_percent.rds")

#only three reps for neuron
metric_percent_neuron <- metric_percent[metric_percent$Rep %in% c(1:3), ]
metric_percent_neuron <- metric_percent_neuron[grepl("Human", metric_percent_neuron$Cluster), ]
metric_percent_neuron$Data <- "Brain \n(neurons)"

# metric_all <- rbind(AMP, 
#                     typist,
#                     tabula,
#                     metric_percent_neuron)
# 
# metric_all$Data <- factor(metric_all$Data,
#                           levels = c("AMP", "Cell typist", "Tabula", "Brain \n(neurons)"))
# 
# metric_all$Method <- factor(metric_all$Method, levels = c("CellGeometry", "MuSiC", "DWLS", "LinDeconSeq"))

metric_allv2 <- rbind(AMP,
                      typist, 
                      tabula,
                      metric_percent[metric_percent$Rep %in% c(1:3),])

metric_allv2$Data[metric_allv2$Data == "Brain"] <- "Human Brain\nCell Atlas"
metric_allv2$Data[metric_allv2$Data == "Tabula"] <- "Tabula Sapiens"

metric_allv2$Data <- factor(metric_allv2$Data,
                            levels = c("AMP", "Cell typist", "Tabula Sapiens", "Human Brain\nCell Atlas"))

metric_allv2$Data_info <- NA
metric_allv2$Data_info[metric_allv2$Data == "AMP" & 
                         metric_allv2$Method == "CellGeometry" &
                         metric_allv2$Cluster == "SC-F1" & metric_allv2$Rep == 1] <- "Synovium\n10,099 cells"
metric_allv2$Data_info[metric_allv2$Data == "Cell typist" & 
                         metric_allv2$Method == "CellGeometry" &
                         metric_allv2$Cluster == "DC" & metric_allv2$Rep == 1] <- "Blood\n27,620 cells"
metric_allv2$Data_info[metric_allv2$Data == "Tabula Sapiens" & 
                         metric_allv2$Method == "CellGeometry" &
                         metric_allv2$Cluster == "fibroblast" & metric_allv2$Rep == 1] <- "75 tissues\n1,136,218 cells"
metric_allv2$Data_info[metric_allv2$Data == "Human Brain\nCell Atlas" & 
                         metric_allv2$Method == "CellGeometry" &
                         metric_allv2$Cluster == "fibroblast" & metric_allv2$Rep == 1] <- "Brain\n3,369,219 cells"

library(ggplot2)
library(ggh4x)

#with removal of outliers
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

library(dplyr)

metric_allv2 <- metric_allv2 %>% group_by(Method, Data) %>%
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

metric_allv2$RMSEv3 <- NA
metric_allv2$RMSEv3[metric_allv2$Data == "AMP" & metric_allv2$RMSE < 6] <- metric_allv2$RMSE[metric_allv2$Data == "AMP" & metric_allv2$RMSE < 6] 
metric_allv2$RMSEv3[metric_allv2$Data == "Cell typist" & metric_allv2$RMSE < 7] <- metric_allv2$RMSE[metric_allv2$Data == "Cell typist" & metric_allv2$RMSE < 7]
metric_allv2$RMSEv3[metric_allv2$Data == "Tabula Sapiens" & metric_allv2$RMSE < 1.14 ] <- metric_allv2$RMSE[metric_allv2$Data == "Tabula Sapiens" & metric_allv2$RMSE < 1.14]
metric_allv2$RMSEv3[metric_allv2$Data == "Human Brain\nCell Atlas" & metric_allv2$RMSE < 1.8 ] <- metric_allv2$RMSE[metric_allv2$Data == "Human Brain\nCell Atlas" & metric_allv2$RMSE < 1.8]

metric_allv2$Method <- factor(metric_allv2$Method,
                              levels = c("CellGeometry",
                                         "MuSiC",
                                         "DWLS",
                                         "LinDeconSeq"))

ggplot(data = metric_allv2) +
  geom_jitter(aes(x = Method, y = RMSEv3, color = Method), alpha = 0.2, 
              width = 0.2)+
  stat_summary(aes(x = Method, y = RMSE), 
               alpha = 0.5,
               fun.data = calc_boxplot_stat, geom="boxplot") + 
  #geom_text(aes(x = Method, y = Inf, label = Data_info),
  #          hjust = 0, nudge_x = -0.5, size = 3, vjust = 1) +
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "LinDeconSeq" = "#00B0F6",
                                "DWLS" = "#E76BF3"),
                     guide = "none") +
  labs(x = "", y = "RMSE") +
  coord_cartesian(clip = "off") +
  facet_grid2(~Data, scales = "free", space = "free_x", independent = "y")+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(1.5, "cm"))


metric_allv2$Rsqv3 <- NA
metric_allv2$Rsqv3[metric_allv2$Data == "AMP" & metric_allv2$Rsq > -2.2] <- metric_allv2$Rsq[metric_allv2$Data == "AMP" & metric_allv2$Rsq > -2.2] 
metric_allv2$Rsqv3[metric_allv2$Data == "Cell typist" & metric_allv2$Rsq > -5.9] <- metric_allv2$Rsq[metric_allv2$Data == "Cell typist" & metric_allv2$Rsq > -5.9] 
metric_allv2$Rsqv3[metric_allv2$Data == "Tabula Sapiens" & metric_allv2$Rsq > -6.2] <- metric_allv2$Rsq[metric_allv2$Data == "Tabula Sapiens" & metric_allv2$Rsq > -6.2] 
metric_allv2$Rsqv3[metric_allv2$Data == "Human Brain\nCell Atlas" & metric_allv2$Rsq > -8.85] <- metric_allv2$Rsq[metric_allv2$Data == "Human Brain\nCell Atlas" & metric_allv2$Rsq > -8.85] 

ggplot(data = metric_allv2) +
  geom_jitter(aes(x = Method, y = Rsqv3, color = Method), alpha = 0.2, 
              width = 0.2)+
  stat_summary(aes(x = Method, y = Rsq), 
               alpha = 0.5,
               fun.data = calc_boxplot_stat, geom="boxplot") + 
  # geom_text(aes(x = Method, y = -Inf, label = Data_info),
  #           hjust = 0, nudge_x = -0.5, size = 3, vjust = -0.5) +
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "LinDeconSeq" = "#00B0F6",
                                "DWLS" = "#E76BF3"),
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
        panel.spacing.x = unit(1.5, "cm"))

