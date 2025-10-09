# Evaluating Tabula deconvolution metrics with similarity 

####RMSE vs similarity####

library(cellGeometry)

metric_percent <- readRDS("tabula_metric_percent.rds")

mk <- readRDS("tabula_markers_dualmeans_rerun.rds")

tabula_mkv2 <- updateMarkers(mk, nsubclass = 200, expfilter = 0.5)

similarity_out <- function(mk, data){
  cs <- cos_similarity(mk)
  
  diag(cs) <- 0
  
  data.frame("Cluster" = names(colMeans(cs)),
             "Mean" = colMeans(cs),
             "Max" = colMaxs(cs),
             "Data" = data)
  
}

cs_df <- similarity_out(tabula_mkv2, "Tabula Sapiens")

library(dplyr)

metric_mean <- metric_percent %>% group_by(Cluster, Method) %>%
  summarise(RMSE_mean = mean(RMSE),
            Rsq_mean = mean(Rsq),
            N = n(),
            RMSE_SD = sd(RMSE),
            Rsq_SD = sd(Rsq),
            RMSE_SEM = RMSE_SD/N,
            Rsq_SEM = Rsq_SD/N)


metric_mean$cs_max <- cs_df$Max[match(metric_mean$Cluster,
                                      cs_df$Cluster)]

library(ggrepel)
library(ggnewscale)
library(ggtext)

max_sim <- function(method, celltypes){
  df <- metric_mean[metric_mean$Method %in% method, ]
  
  df$label <- NA
  df$label[df$Cluster %in% celltypes] <- as.vector(df$Cluster[df$Cluster %in% celltypes])
  
  if(length(method) == 1){
    nls_fit <- nls(RMSE_mean ~ a + c / (b - cs_max), start = list(a = 0.1, b = 1, c = 1),
                   data = df)
    print(coef(nls_fit))
    rsq <- modelr::rsquare(nls_fit, df)
    print(rsq)
    rsq_lab <- paste0("R<sup>2</sup> = ", signif(rsq, 2))
  }
  else{rsq_lab <- NULL}
  
  ggplot(df, aes(x = cs_max, y = RMSE_mean, color = Method, group = Method)) +
    geom_point(alpha = 0.55, size = 1.5) +
    scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                  "MuSiC" = "#00BF7D",
                                  "LinDeconSeq" = "#00B0F6",
                                  "DWLS" = "#E76BF3"),
                       guide = "none") +
    new_scale_color()+
    geom_smooth(aes(color = Method),
                method = "nls",
                formula = y ~ a + c / (b - x),
                se = FALSE, 
                method.args = list(start = list(a = 0.1, b = 1, c = 1)), 
                linewidth = 0.5,  fullrange = TRUE) +
    scale_x_continuous(limits = c(min(df$cs_max) -0.02, max(df$cs_max) + 0.0015), expand = c(0, 0))+
    scale_color_manual(values = c("CellGeometry" = "#b51207", 
                                  "MuSiC" = "#05754e"),
                       guide = "none") +
    coord_cartesian(clip = "off") +
    labs(x = "Max similarity",
         y = "RMSE mean",
         title = method,
         subtitle = rsq_lab) +
    theme_classic() +
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(size = 11),
          plot.subtitle = element_markdown(),
          plot.margin = margin(1.5, 10, 1.5, 1.5))
}

metric_mean$Cluster <- as.vector(metric_mean$Cluster)

metric_mean$Cluster[metric_mean$Cluster == "CD4-positive, alpha-beta T cell"] <- "CD4 T cell"
metric_mean$Cluster[metric_mean$Cluster == "naive thymus-derived CD4-positive, alpha-beta T cell"] <- "naive thymus-derived CD4 T cell"

cellgeo_tab <- max_sim("CellGeometry", 
                       c("CD4 T cell",
                         "naive thymus-derived CD4 T cell",
                         "fibroblast",
                         "monocyte",
                         "pancreatic acinar cell",
                         "oocyte", "platelet", "erythrocyte")) 

music_tab <- max_sim("MuSiC",
                     c("macrophage", "basal cell"))

pdf("RMSE_cs_max_tabula_updateMarkers_margins_fittedline_nolabs_correct.pdf", width = 11.75, height = 3.96)
ggarrange(cellgeo_tab, music_tab, comb_tab, align = "hv", ncol = 3)
dev.off()

####specificity plot####

source("specificity_plot.R")

pdf("tabula_specificity_updateMarkersv2.pdf",
    width =  13.5, height = 7.3)
ggarrange(specificity_plotv2(tabula_mkv2, "oocyte",
                             add_labels = c("GDF9", "KPNA7", "PADI6",
                                            "TUBB8", "ZP3"), force = 5,
                             min.segment.length = unit(0, 'lines')),
          specificity_plotv2(tabula_mkv2, "platelet",
                             add_labels = c("ITGB3", "TUBB1", "GP9", "PF4V1",
                                            "PPBP"),
                             force = 5,
                             min.segment.length = unit(0, 'lines')),
          specificity_plotv2(tabula_mkv2, "pancreatic acinar cell", 
                             add_labels = c("PRSS1", "CPA1", "PNLIP",
                                            "CPB1", "CTRB1"),
                             force = 8,
                             min.segment.length = unit(0, 'lines')),
          specificity_plotv2(tabula_mkv2, "erythrocyte",
                             add_labels = c("HBA2", "HBA1", "ALAS2", "HBM",
                                            "HBQ1")),
          specificity_plotv2(tabula_mkv2, "CD4-positive, alpha-beta T cell",
                             add_labels = c("ZNF831", "LTB", "ICOS"),
                             force = 5,
                             min.segment.length = unit(0, 'lines'),
                             max.overlaps = 40) + ggtitle("CD4 T cell"),
          specificity_plotv2(tabula_mkv2, "naive thymus-derived CD4-positive, alpha-beta T cell",
                             add_labels = c("LEF1", "LTB", "TCF7", "CCR7", "CD28"),
                             force = 5,
                             min.segment.length = unit(0, 'lines'),
                             max.overlaps = 40) + ggtitle("naive thymus-derived CD4 T cell"),
          specificity_plotv2(tabula_mkv2, "fibroblast",  
                             add_labels = c("TWIST2", "DPT"), box.padding = 0,
                             force_pull = 0),
          specificity_plotv2(tabula_mkv2, "thymic fibroblast type 2",
                             add_labels = c("VEGFD", "PI16", "RAMP2-AS1"),
                             force = 5,
                             min.segment.length = unit(0, 'lines'),
                             max.overlaps = 40), align = "hv",
          ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()
