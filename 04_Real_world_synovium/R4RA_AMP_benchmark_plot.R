# comparing the deconvolution methods for R4RA synovium (AMP dataset reference)

r4ra_ngroup <- readRDS("Cellgeo_AMP_rawcounts_countT_eq200_ngroup.rds")

cellgeo <- r4ra_ngroup$`Ngroup = 10`
rm(r4ra_ngroup)


clindata <- readRDS("/media/gcpeac/Myles/Deconvolution/metadata.R4RA.rds")
col_order <- c(paste0('SC-F', 1:4), paste0('SC-M', 1:4), paste0('SC-T', 1:6),
               paste0('SC-B', 1:4))

music2fit <- readRDS("music2_AMP_rawcounts.rds")
load("DWLS_AMP.rdata")
load("LinDeconSeq_AMPv2_rawcounts.rdata")

cellgeo_percent <- as.data.frame(cellgeo$nest_percent)
music_percent <- as.data.frame(music2fit$Est.prop.weighted * 100)
DWLS_percent <- as.data.frame(t(allCounts_DWLS  * 100))
Lin_percent <- as.data.frame(fractionsv2 * 100)
colnames(Lin_percent) <- gsub("\\.", "-", colnames(Lin_percent))
rownames(Lin_percent) <- gsub("\\.", "-", rownames(Lin_percent))

library(tidyr)
load("/media/gcpeac/R4RA/DB_releases/R4RA_270821.RData")
r4ra.hist <- R4RA_db_v202108$lab.db
r4ra.histBL <- r4ra.hist[r4ra.hist$Visit == 2, ]
r4ra.histv7 <- r4ra.hist[r4ra.hist$Visit == 7, ]

arrange_data <- function(df, values, method){
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = values))
  test$Method <- method
  test$PathotypeBL <- as.vector(clindata$Pathotype.V2[match(test$Patient,
                                                            clindata$Seq_ID.V2)])
  test$PathotypeBL <- factor(test$PathotypeBL,
                             levels = c("Ungraded", "Fibroid", "Myeloid", "Lymphoid"))
  test$Pathotype <- as.vector(r4ra.meta$Pathotype[match(test$Patient,
                                                        r4ra.meta$noPrefix)])
  test$Pathotype[test$Pathotype == "TBC"] <- "Ungraded"
  test$Pathotype <- factor(test$Pathotype,
                           levels = c("Ungraded", "Fibroid", "Myeloid", "Lymphoid"))
  
  test$Cluster <- factor(test$Cluster)
  test$Celltype <- NA
  test$Celltype[grepl("F", test$Cluster)] <- "Fibroblast"
  test$Celltype[grepl("M", test$Cluster)] <- "Macrophage"
  test$Celltype[grepl("B", test$Cluster)] <- "B cell"
  test$Celltype[grepl("T", test$Cluster)] <- "T cell"
  
  test$Celltype <- factor(test$Celltype,
                          levels = c("Fibroblast", "Macrophage",
                                     "T cell", "B cell"))
  test$ID <- r4ra.meta$Patient.I.D.[match(test$Patient,
                                          r4ra.meta$noPrefix)]
  
  test$CD138 <- NA
  test$CD138[grepl("-Baseline", test$Patient)] <- r4ra.histBL$CD138[match(test$ID[grepl("-Baseline", test$Patient)],
                                                                          r4ra.histBL$Patient.I.D.)]
  
  test$CD138[!grepl("-Baseline", test$Patient)] <- r4ra.histv7$CD138[match(test$ID[!grepl("-Baseline", test$Patient)],
                                                                           r4ra.histv7$Patient.I.D.)]
  
  test$CD3 <- NA
  test$CD3[grepl("-Baseline", test$Patient)] <- r4ra.histBL$CD3[match(test$ID[grepl("-Baseline", test$Patient)],
                                                                      r4ra.histBL$Patient.I.D.)]
  
  test$CD3[!grepl("-Baseline", test$Patient)] <- r4ra.histv7$CD3[match(test$ID[!grepl("-Baseline", test$Patient)],
                                                                       r4ra.histv7$Patient.I.D.)]
  
  test$CD20 <- NA
  test$CD20[grepl("-Baseline", test$Patient)] <- r4ra.histBL$CD20[match(test$ID[grepl("-Baseline", test$Patient)],
                                                                        r4ra.histBL$Patient.I.D.)]
  
  test$CD20[!grepl("-Baseline", test$Patient)] <- r4ra.histv7$CD20[match(test$ID[!grepl("-Baseline", test$Patient)],
                                                                         r4ra.histv7$Patient.I.D.)]
  
  test$CD68SL <- NA
  test$CD68SL[grepl("-Baseline", test$Patient)] <- r4ra.histBL$CD68SL[match(test$ID[grepl("-Baseline", test$Patient)],
                                                                            r4ra.histBL$Patient.I.D.)]
  
  test$CD68SL[!grepl("-Baseline", test$Patient)] <- r4ra.histv7$CD68SL[match(test$ID[!grepl("-Baseline", test$Patient)],
                                                                             r4ra.histv7$Patient.I.D.)]
  
  test$CD68L <- NA
  test$CD68L[grepl("-Baseline", test$Patient)] <- r4ra.histBL$CD68L[match(test$ID[grepl("-Baseline", test$Patient)],
                                                                          r4ra.histBL$Patient.I.D.)]
  
  test$CD68L[!grepl("-Baseline", test$Patient)] <- r4ra.histv7$CD68L[match(test$ID[!grepl("-Baseline", test$Patient)],
                                                                           r4ra.histv7$Patient.I.D.)]
  
  test
}

comb <- rbind(arrange_data(cellgeo_percent, "Percentage", "CellGeometry"),
              arrange_data(music_percent, "Percentage", "MuSiC"),
              arrange_data(DWLS_percent, "Percentage", "DWLS"),
              arrange_data(Lin_percent, "Percentage", "LinDeconSeq"))

comb$Cluster <- factor(comb$Cluster,
                       levels = col_order)

comb$Method <- factor(comb$Method,
                      levels = c("CellGeometry",
                                 "MuSiC",
                                 "DWLS",
                                 "LinDeconSeq"))

cellgeo_group_output <- arrange_data(as.data.frame(cellgeo$group$output), "Output", "CellGeometry")
cellgeo_nest_output <- arrange_data(as.data.frame(cellgeo$nest_output), "Output", "CellGeometry")

####CD markers####

combv2 <- comb
combv2$Percentage[combv2$Method == "DWLS" & combv2$Percentage == 0] <- NA

CD_plot <- function(df, cell, CD, Method, param = "Percentage"){
  test <- df
  test$CD <- test[ , CD]
  test <- test[is.na(test$CD) == FALSE, ]
  test <- test[test$Cluster == cell, ]
  test$Param <- test[ , param]
  
  stats <- cor.test(test$CD, test$Param, method = "spearman")
  
  if(param == "Output"){
    ylab <- cell
  }
  else{
    ylab <- paste(cell, "(%)")
  }
  
  ggplot(test, aes(x = CD, y = Param)) +
    geom_point(alpha = 0.8) +
    stat_smooth(method = "lm", col = "red") +
    labs(y = ylab, x = CD,
         title = Method, 
         subtitle = paste0("r = ", formatC(stats$estimate,2),
                           " P = ", formatC(stats$p.value,2))) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 12))
  
}

#cell groups of interest

CD_group_plot <- function(df, cell, CD, Method, param = "Percentage"){
  test <- df
  test$CD <- test[ , CD]
  test <- test[is.na(test$CD) == FALSE, ]
  test <- test[test$Celltype == cell, ]
  test$Param <- test[ , param]
  
  testv2 <- as.data.frame(test[ , c("Patient", "Celltype", "CD", "Param")] %>%
                            group_by(Patient) %>%
                            mutate(Sum = sum(Param)))
  
  testv2 <- distinct(testv2[ , c("Patient", "Celltype", "CD", "Sum")])
  stats <- cor.test(testv2$CD, testv2$Sum, method = "spearman")
  
  if(param == "Output"){
    ylab <- cell
  }
  else{
    ylab <- paste(cell, "(%)")
  }
  
  ggplot(testv2, aes(x = CD, y = Sum)) +
    geom_point(alpha = 0.8) +
    stat_smooth(method = "lm", col = "red") +
    labs(y = ylab,
         x = CD, 
         title = Method, 
         subtitle = paste0("r = ", formatC(stats$estimate,2),
                           " P = ", formatC(stats$p.value,2))) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 12))
}

cellgeoB <- list()
cellgeoB[["1"]] <- CD_plot(cellgeo_group_output, "B cell" , "CD20", "CellGeometry", "Output")

plotlistB <- lapply(c("MuSiC", "DWLS", "LinDeconSeq"), function(x){
  CD_group_plot(combv2[combv2$Method == x, ], "B cell", "CD20", x, "Percentage")
})

plotlistB <- c(cellgeoB,
               plotlistB)

cellgeoM <- list()
cellgeoM[["1"]] <- CD_plot(cellgeo_group_output, "Monocyte" , "CD68SL", "CellGeometry", "Output") +
  ylab("Macrophage")

plotlistM <- lapply(c("MuSiC", "DWLS", "LinDeconSeq"), function(x){
  CD_group_plot(combv2[combv2$Method == x, ], "Macrophage", "CD68SL", x, "Percentage")
})

plotlistM <- c(cellgeoM,
               plotlistM)

cellgeoT <- list()
cellgeoT[["1"]] <- CD_plot(cellgeo_group_output, "T cell" , "CD3", "CellGeometry", "Output") 

plotlistT <- lapply(c("MuSiC", "DWLS", "LinDeconSeq"), function(x){
  CD_group_plot(combv2[combv2$Method == x, ], "T cell", "CD3", x, "Percentage")
})

plotlistT <- c(cellgeoT,
               plotlistT)

cellgeoP <- list()
cellgeoP[["1"]] <- CD_plot(cellgeo_nest_output, "SC-B4" , "CD138", "CellGeometry", "Output")

plotlistP <- lapply(c("MuSiC", "DWLS", "LinDeconSeq"), function(x){
  CD_plot(combv2[combv2$Method == x, ], "SC-B4", "CD138", x)
})

plotlistP <- c(cellgeoP,
               plotlistP)

plotlist_comb <- c(plotlistM,
                   plotlistT,
                   plotlistB,
                   plotlistP)

ggarrange(plotlist = plotlist_comb, ncol = 4, nrow = 4,
          align = "v")

####stacked barplots####

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

order_method <- "CA"

#for CA need to remove constant NA/zero

NA_rows <- rowSums(DWLS_percent) == 0

DWLS_percentv2 <- DWLS_percent[!NA_rows, ]

NA_cols <- colSums(Lin_percent) == 0

Lin_percentv2 <- Lin_percent[ , !NA_cols]

patho_col <- c("Lymphoid" = "steelblue",
               "Myeloid" = "firebrick",
               "Fibroid" = "palegreen",
               "Ungraded" = "grey80")

col_scheme <- c(amp_scheme, patho_col)

bar_outv4 <- function(percent, method, order_method, rev = FALSE){
  #Requires a nonnegative matrix.
  
  o <- seriate(percent, method = order_method)
  
  comb_sub <- comb[comb$Method == method, ]
  
  patient <- rownames(percent)[o[[1]]]
  
  if(rev == TRUE){
    patient <- rev(patient)
  }
  
  comb_sub$Patient <- factor(comb_sub$Patient,
                             levels = patient)
  
  comb_sub <- comb_sub[comb_sub$Percentage > 0, ]
  
  df <- data.frame("Patient" = unique(comb_sub$Patient))
  df$Patient <- factor(df$Patient, levels = patient)
  df$Pathotype <- comb_sub$Pathotype[match(df$Patient,
                                           comb_sub$Patient)]
  df$Total = 1.13 
  
  ggplot(comb_sub) +
    geom_col(aes(x = Patient, y = Percentage, fill = Cluster), 
             position = "fill", width = 1) +
    geom_tile(data = df, aes(x = Patient, y = Total, fill = Pathotype, color = Pathotype),
              height = 0.11, width = 1)+ #was 0.05 originally
    scale_fill_manual(values = col_scheme, breaks = names(amp_scheme)) +
    scale_color_manual(values = col_scheme, breaks = names(amp_scheme), guide = "none") +
    scale_y_continuous(labels = scales::percent_format(suffix = "")) +
    coord_cartesian(clip = "off", ylim = c(0, 1)) +
    labs(x = "", y = "Cell (%)", title = method)+
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 8),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 8),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 10, vjust = 4.5)) #+
  #guides(fill = guide_legend(ncol = 9))
}

fig <- ggarrange(bar_outv4(cellgeo_percent, "CellGeometry", order_method, rev = TRUE) +
                   theme(legend.position = "none"),
                 bar_outv4(music_percent, "MuSiC", order_method, rev = TRUE)+
                   theme(legend.position = "none"),
                 bar_outv4(DWLS_percentv2, "DWLS", order_method, rev = TRUE)+
                   theme(legend.position = "none"),
                 bar_outv4(Lin_percentv2, "LinDeconSeq", order_method)+
                   theme(legend.position = "none"),
                 ncol = 1, nrow = 4)

annotate_figure(fig, top = text_grob("R4RA synovium", size = 10))

####stacked barplots by subclass####

cellgeo_output_mean <- cellgeo_nest_output %>% group_by(Cluster) %>%
  summarise("Mean" = mean(Output, na.rm = TRUE))

cellgeo_output_mean_patho <- as.data.frame(cellgeo_nest_output %>% group_by(Cluster, Pathotype) %>%
                                             summarise("Mean" = mean(Output, na.rm = TRUE)))

cellgeo_output_mean_patho <- cellgeo_output_mean_patho[cellgeo_output_mean_patho$Pathotype != "Ungraded", ]
cellgeo_output_mean_patho$Celltype <- cellgeo_nest_output$Celltype[match(cellgeo_output_mean_patho$Cluster,
                                                                         cellgeo_nest_output$Cluster)]

comb_mean_patho <- as.data.frame(combv2 %>% group_by(Cluster, Pathotype, Method) %>%
                                   summarise("Mean" = mean(Percentage, na.rm = TRUE)))

comb_mean_patho <- comb_mean_patho[comb_mean_patho$Pathotype != "Ungraded", ]
comb_mean_patho$Celltype <- combv2$Celltype[match(comb_mean_patho$Cluster,
                                                  combv2$Cluster)]

comb_mean_patho$Pathotype <- factor(comb_mean_patho$Pathotype,
                                    levels = c("Fibroid", "Myeloid", "Lymphoid"))

patho_stackv2 <- function(method, value = "Percentage"){
  if(value == "Output"){
    sub <- cellgeo_output_mean_patho
  } else{
    sub <- comb_mean_patho[comb_mean_patho$Method == method, ]
  }
  
  if(value == "Output"){
    ylab = "Cell freq"
  } else{
    ylab = "Cell (%)"
  }
  
  ggplot() + 
    geom_col(data = sub, aes(x = Pathotype, y = Mean, fill = Cluster), 
             width = 0.5, color = "black") + 
    ggtitle(method)+
    labs(x = "", y = ylab) +
    facet_wrap(~Celltype, 
               #nrow = 4,
               nrow = 1, ncol = 4,
               scales = "free_y")+
    scale_fill_manual(values = amp_scheme,
                      breaks = names(amp_scheme),
                      name = "",
                      guide = "none")+
    theme_classic()+
    theme(axis.text = element_text(color = "black", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12),
          strip.text = element_blank(),
          plot.title = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          strip.background = element_blank())
}

ggarrange(
  patho_stackv2("CellGeometry", "Output"),
  patho_stackv2("MuSiC") ,
  patho_stackv2("DWLS") + ylab(""),
  patho_stackv2("LinDeconSeq") + ylab (""),
  ncol = 1, nrow = 4)

####Table####

#detailing the mean cell proportions and its range 

library(dplyr)

stats_df <- combv2 %>% group_by(Cluster, Method) %>%
  summarise("Mean" = mean(Percentage, na.rm = TRUE),
            "Min" = min(Percentage, na.rm = TRUE),
            "Max" = max(Percentage, na.rm = TRUE))

stats_df$Mean <- format(round(stats_df$Mean, 2), nsmall = 2)
stats_df$Min <- format(round(stats_df$Min, 2), nsmall = 2)
stats_df$Max <- format(round(stats_df$Max, 2), nsmall = 2)

stats_df$Range <- paste0(trimws(stats_df$Mean),
                         " (",
                         trimws(stats_df$Min),
                         " - ",
                         trimws(stats_df$Max),
                         ")")

stats_df_wide <- pivot_wider(stats_df,
                             id_cols = "Cluster",
                             names_from = "Method",
                             values_from = "Range")

# write.table(stats_df_wide, file = "R4RA_prop_df.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)








