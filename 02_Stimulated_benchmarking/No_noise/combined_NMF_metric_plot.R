####Background####

# plotting metric plots for NMF and cellGeometry

####Organise data####

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

arrange_metricv2 <- function(list, times){
  df <- data.frame()
  for(i in names(list)){
    for(x in names(list[[i]])){
      temp <- as.data.frame(list[[i]][[x]]$metrics_percent)
      temp <- add_column(temp,
                         "Cluster" = rownames(temp),
                         .before = 1)
      temp$Times <- times
      temp$Rep <- as.numeric(gsub("Rep ", "", x))
      temp$Method <- i
      df <- rbind(df,
                  temp)
    }
  }
  df
}

#####AMP#####

AMP_cellgeo <- readRDS("../../../AMP/Dirichlet/cellgeo_dirichlet_output500.rds")
AMP_NMF <- readRDS("../../../AMP/Dirichlet/NMF_dirichlet_output.rds")
AMP_NMFalt <- readRDS("../../../AMP/Dirichlet/NMF_dirichlet_output_alt.rds")
AMP_rcppML <- readRDS("../../../AMP/Dirichlet/RcppML_dirichlet_output.rds")

AMP_metric <- rbind(arrange_metric(AMP_cellgeo, "CellGeometry"),
                    arrange_metric(AMP_NMF, "Brunet"),
                    arrange_metricv2(AMP_NMFalt, 30),
                    arrange_metric(AMP_rcppML, "RcppML"))

#####typist####

typist_cellgeo <- readRDS("../../../Cell_typist/workstation2/Dirichlet/cellgeo_dirichlet_output.rds")
typist_NMF <- readRDS("../../../Cell_typist/workstation2/Dirichlet/NMF_dirichlet_output.rds")
typist_NMFalt <- readRDS("../../../Cell_typist/workstation2/Dirichlet/NMF_dirichlet_output_alt.rds")
typist_rcppML <- readRDS("../../../Cell_typist/workstation2/Dirichlet/RcppML_dirichlet_output.rds")

typist_metric <- rbind(arrange_metric(typist_cellgeo, "CellGeometry"),
                      arrange_metric(typist_NMF, "Brunet"),
                      arrange_metricv2(typist_NMFalt, 30),
                      arrange_metric(typist_rcppML, "RcppML"))

#####tabula#####

tabula_cellgeo <- readRDS("../../../Tabula/workstation2/cellgeo_dirichlet_output_500.rds")
tabula_NMF <- readRDS("../../../Tabula/workstation2/NMF_dirichlet_output.rds")
tabula_NMFalt <- readRDS("../../../Tabula/workstation2/NMF_dirichlet_output_alt.rds")
tabula_rcppML <- readRDS("../../../Tabula/workstation2/RcppML_dirichlet_output.rds")

tabula_metric <- rbind(arrange_metric(tabula_cellgeo, "CellGeometry"),
                       arrange_metric(tabula_NMF, "Brunet"),
                       arrange_metricv2(tabula_NMFalt, 3),
                       arrange_metric(tabula_rcppML, "RcppML"))

rm(list = setdiff(ls(), c("AMP_metric", "typist_metric",
                          "tabula_metric", "arrange_metric",
                          "arrange_metricv2")))

#####brain#####

brain_cellgeo <- readRDS("../../cellgeo_dirichlet_output.rds")
brain_NMF <- readRDS("../../NMF_dirichlet_output.rds")
brain_NMFalt <- readRDS("../../NMF_dirichlet_output_alt.rds")
brain_rcppML <- readRDS("../../RcppML_dirichlet_output.rds")

brain_metric <- rbind(arrange_metric(brain_cellgeo, "CellGeometry"),
                       arrange_metric(brain_NMF, "Brunet"),
                       arrange_metricv2(brain_NMFalt, 1),
                       arrange_metric(brain_rcppML, "RcppML"))

rm(list = setdiff(ls(), c("AMP_metric", "typist_metric",
                          "tabula_metric","brain_metric",
                          "arrange_metric",
                          "arrange_metricv2")))

####combine####

AMP_metric$Data <- "AMP"
typist_metric$Data <- "Cell typist"
tabula_metric$Data <- "Tabula Sapiens"
brain_metric$Data <- "Human Brain\nCell Atlas"

metric_all <- rbind(AMP_metric,
                    typist_metric, 
                    tabula_metric,
                    brain_metric[brain_metric$Rep %in% c(1:3),])

metric_all$Method[metric_all$Method == "lee"] <- "Lee"

metric_all$Method <- factor(metric_all$Method,
                            levels = c("CellGeometry",
                                       "Brunet",
                                       "KL",
                                       "Lee",
                                       "RcppML"))

metric_all$Data <- factor(metric_all$Data,
                          levels = c("AMP", "Cell typist", "Tabula Sapiens", "Human Brain\nCell Atlas"))

metric_all$Data_info <- NA
metric_all$Data_info[metric_all$Data == "AMP" & 
                       metric_all$Method == "CellGeometry" &
                       metric_all$Cluster == "SC-F1" & metric_all$Rep == 1] <- "Synovium\n10,099 cells"
metric_all$Data_info[metric_all$Data == "Cell typist" & 
                       metric_all$Method == "CellGeometry" &
                       metric_all$Cluster == "DC" & metric_all$Rep == 1] <- "Blood\n27,620 cells"
metric_all$Data_info[metric_all$Data == "Tabula Sapiens" & 
                       metric_all$Method == "CellGeometry" &
                       metric_all$Cluster == "fibroblast" & metric_all$Rep == 1] <- "75 tissues\n1,136,218 cells"
metric_all$Data_info[metric_all$Data == "Human Brain\nCell Atlas" & 
                       metric_all$Method == "CellGeometry" &
                       metric_all$Cluster == "fibroblast" & metric_all$Rep == 1] <- "Brain\n3,369,219 cells"


library(ggh4x)
library(ggplot2)

ggplot(data = metric_all) +
  geom_jitter(aes(x = Method, y = RMSE, color = Method), alpha = 0.2, 
              width = 0.2)+
  geom_boxplot(aes(x = Method, y = RMSE), 
               alpha = 0.5, outlier.shape = NA) + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "Brunet" = "grey60",
                                "KL" = "burlywood3",
                                "Lee" = "slategray3",
                                "RcppML" = "thistle3"),
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
        panel.spacing.x = unit(1.5, "cm"))


ggplot(data = metric_all) +
  geom_jitter(aes(x = Method, y = Rsq, color = Method), alpha = 0.2, 
              width = 0.2)+
  geom_boxplot(aes(x = Method, y = Rsq), 
               alpha = 0.5, outlier.shape = NA) + 
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "Brunet" = "grey60",
                                "KL" = "burlywood3",
                                "Lee" = "slategray3",
                                "RcppML" = "thistle3"),
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

