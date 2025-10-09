# plots to compare deconvolution of PEAC blood deconvolution 

path <- "/media/gcpeac/Rachel/Packages/cellGeometry_paper/PEAC/Typist_blood/"

music <- readRDS(paste0(path, "music2_PEAC_typist_rawcounts.rds"))
load(paste0(path, "DWLS_PEAC_typist_rawcounts.rdata"))
load(paste0(path, "LinDeconSeq_PEAC_typist.rdata"))

library(cellGeometry)

cellgeo <- readRDS("Cellgeo_PEAC_merge.rds")

mk <- readRDS("merge_markers_hd5_collapse.rds")

col_order <- c(
  "Memory B cells", "Naive B cells",
  "Plasma cells",
  "DC", "DC2", "pDC", 
  "Classical monocytes", "Non-classical monocytes",
  "Macrophages",
  "Follicular helper T cells", "MAIT cells", "Regulatory T cells",
  "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells",
  "Tem/Effector cytotoxic T cells", "Tem/Effector helper T cells",
  "Cytotoxic T cells", "Helper T cells", "Gamma delta T cells",
  "Cycling T cells",
  "CD16+ NK cells", "CD16- NK cells", "NK cells", "ILC",
  #"Erythrocyte",
  #"Platelet",
  "Mast cells",
  "Basophil",
  "Neutrophil",
  "Early MK", "MEMP", "Promyelocytes",
  "HSC/MPP", 
  "GMP")

col_scheme <- c("Memory B" = "skyblue1",
                "Naive B" = "lightblue",
                "Plasma" = "royalblue4",
                "DC" = "sienna4",
                "DC2" = "sienna3",
                "pDC" = "sienna1",
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
                "Cytotoxic T" = "plum4", 
                "Helper T" = "purple3", 
                "Gamma delta T" = "orchid",
                "Cycling T" = "#490092",
                "CD16+ NK" = "aquamarine",
                "CD16- NK" = "aquamarine3",
                "NK" = "aquamarine4",
                "ILC" = "lightseagreen",
                #"Erythrocyte" = "indianred",
                #"Platelet" = "indianred4",
                "Mast" = "orange",
                "Basophil" = "darkgoldenrod3",
                "Neutrophil" = "gold",
                "Early MK" = "gray75",
                "MEMP" = "gray90", 
                "Promyelocytes" = "gray30",
                "HSC/MPP" = "gray50",
                "GMP" = "gray20")

cell_groups <- as.data.frame(mk$cell_table)
colnames(cell_groups) <- "Group"
cell_groups$Group <- as.vector(cell_groups$Group)
cell_groups$Group[cell_groups$Group == "granulocyte"] <- "Granulocyte"

####organise data####

prop_arrange <- function(df){
  prop <- df$Est.prop.weighted * 100
  #common <- intersect(colnames(cellgeo$subclass$percent), colnames(prop))
  
  #prop <- prop[ , colnames(prop) %in% common]
  
  missing <- setdiff(colnames(cellgeo$subclass$percent), colnames(prop))
  
  propdf <- as.data.frame(prop)
  
  if(length(missing) > 0){
    for(i in missing){
      propdf[ , i] <- NA
    }
  }
  
  #propdf[ , match(colnames(cellgeo$subclass$percent), colnames(propdf))]
  propdf
}

cellgeo_percent <- as.data.frame(cellgeo$subclass$percent)

temp <- music$Est.prop.weighted

colnames(temp)[colSums(temp) == 0]
# [1] "Regulatory T cells"          "Memory B cells"              "Tcm/Naive helper T cells"    "Tcm/Naive cytotoxic T cells" "CD16- NK cells"             
# [6] "Naive B cells"               "pDC"                         "DC"                          "Helper T cells"              "Mast cells"                 
# [11] "ILC"                         "MEMP"                        "HSC/MPP"        

setdiff(colnames(cellgeo_percent), colnames(temp))
# [1] "Follicular helper T cells"      "GMP"                            "Tem/Effector cytotoxic T cells" "basophil"                      
# [5] "neutrophil"  

music_percent <- as.data.frame(prop_arrange(music))

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
  
  print(setdiff(colnames(df), colnames(cellgeo$subclass$percent)))
  df
  #common <- intersect(colnames(cellgeo$subclass$percent), colnames(df))
  #df <- df[ , colnames(df) %in% common]
  #df[ , match(colnames(cellgeo$subclass$percent), colnames(df))]
}

DWLS_percent <- as.data.frame(prop_arrangev2(t(allCounts_DWLS * 100)))
#"Cytotoxic T cells"   "Helper T cells"      "MEMP"                "Promyelocytes"       "gamma delta T cells"

Lin_percent <- as.data.frame(prop_arrangev2(fractions * 100))
#"Cytotoxic T cells"   "Helper T cells"      "MEMP"                "Promyelocytes"       "gamma delta T cells"

colnames(Lin_percent)[colSums(Lin_percent) == 0]
#[1] "CD16+ NK cells"                 "NK cells"                       "DC"                            
# [4] "Cycling T cells"                "Macrophages"                    "HSC/MPP"                       
# [7] "Tem/Effector cytotoxic T cells" "Follicular helper T cells"      "gamma delta T cells" 

library(tidyr)
arrange_data <- function(df, value, method){
  df$Patient <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Patient,
                                     names_to = "Cluster",
                                     values_to = value))
  test$Method <- method
  test$Group <- cell_groups$Group[match(test$Cluster,
                                        rownames(cell_groups))]
  filter <- c("erythrocyte", "platelet", "basophil", "neutrophil")
  test$Cluster[test$Cluster %in% filter] <- stringr::str_to_sentence(test$Cluster[test$Cluster %in% filter])
  #test$Cluster <- factor(test$Cluster, levels = col_order)
  test
}

comb <- rbind(arrange_data(cellgeo_percent, "Percentage", "CellGeometry"),
              arrange_data(music_percent, "Percentage","MuSiC"),
              arrange_data(DWLS_percent, "Percentage","DWLS"),
              arrange_data(Lin_percent, "Percentage","LinDeconSeq"))

comb$Cluster[comb$Cluster == "gamma delta T cells"] <- "Gamma delta T cells"

comb$Cluster <- factor(comb$Cluster,
                       levels = col_order)

comb$Method <- factor(comb$Method,
                      levels = c("CellGeometry",
                                 "MuSiC",
                                 "DWLS",
                                 "LinDeconSeq"))

####stacked barplot####

bar_outv4 <- function(percent, method, order_method){
  #Requires a nonnegative matrix.
  percent[is.na(percent) == TRUE] <- 0
  o <- seriate(percent, method = order_method)
  
  comb_sub <- comb[comb$Method == method, ]
  
  patient <- rownames(percent)[o[[1]]]
  
  comb_sub$Patient <- factor(comb_sub$Patient,
                             levels = patient)
  
  comb_sub$Percentage[is.na(comb_sub$Percentage) == TRUE] <- 0
  comb_sub <- comb_sub[comb_sub$Percentage > 0, ]
  
  comb_sub$Cluster <- as.vector(comb_sub$Cluster)
  comb_sub$Cluster <- gsub(" cells", "", comb_sub$Cluster)
  comb_sub$Cluster <- factor(comb_sub$Cluster,
                             levels = names(col_scheme))
  
  ggplot(comb_sub) +
    geom_bar(aes(x = Patient, y = Percentage, fill = Cluster), 
             position = "fill", stat = "identity", width = 1) +
    scale_fill_manual(values = col_scheme, breaks = names(col_scheme), drop = FALSE)+
    #scale_color_manual(values = col_scheme, breaks = names(col_scheme), drop = FALSE)+
    scale_y_continuous(labels = scales::percent_format(suffix = "")) +
    labs(x = "", y = "Cell (%)", title = method)+
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 8),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 8),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(fill=guide_legend(ncol=2)) 
}

library(seriation)

order_method <- "CA"

NA_cols <- colSums(Lin_percent) == 0
Lin_percentv2 <- Lin_percent[ , !NA_cols]

NA_colsv2 <- colSums(music_percent) == 0 | colSums(is.na(music_percent)) == nrow(music_percent)
music_percentv2 <- music_percent[ , !NA_colsv2]

library(ggplot2)
library(ggpubr)

fig <- ggarrange(bar_outv4(cellgeo_percent, "CellGeometry", order_method) +
                   theme(legend.position = "none"),
                 bar_outv4(music_percentv2, "MuSiC", order_method)+
                   theme(legend.position = "none"),
                 bar_outv4(DWLS_percent, "DWLS", order_method)+
                   theme(legend.position = "none"),
                 bar_outv4(Lin_percentv2, "LinDeconSeq", order_method)+
                   theme(legend.position = "none"),
                 ncol = 2, nrow = 2)

pdf("Basophil_neutro/PEAC_percent_methods.PDF", width = 6.2, height = 3)
annotate_figure(fig, top = text_grob("PEAC blood", size = 10))
dev.off()

####Pie chart####

col_scheme2 <- c("B cells" = "steelblue",
                 "Granulocyte" = "gold",
                 "ILC" = "aquamarine",
                 "Monocytes" = "firebrick1",
                 "Progenitor/misc" = "gray",
                 "T cells" = "purple")

comb$Group[comb$Cluster %in% c("Cytotoxic T cells",
                               "Helper T cells", "Gamma delta T cells")] <- "T cells"

comb$Group[comb$Cluster %in% c("MEMP", "Promyelocytes")] <- "Progenitor/misc"

comb_sum <- comb[ , c("Percentage", "Patient", "Method", "Group")] %>%
  group_by(Group, Patient, Method) %>%
  summarise("Sum" = sum(Percentage, na.rm = TRUE))

comb_meanv2 <- comb_sum %>%
  group_by(Group, Method) %>%
  summarise("Mean" = mean(Sum, na.rm = TRUE))

comb_meanv2$Data <- c("Cell typist")
comb_meanv2$Data[comb_meanv2$Method == "CellGeometry"] <- "Cell typist + Tabula"
comb_meanv2$Data <- factor(comb_meanv2$Data,
                           levels = c("Cell typist + Tabula",
                                      "Cell typist"))

comb_meanv2$Group <- factor(comb_meanv2$Group,
                            levels = c("B cells", "Monocytes",
                                       "T cells", "ILC", 
                                       "Granulocyte",
                                       "Progenitor/misc"))

library(ggh4x)

P <- ggplot()+
  geom_col(data = comb_meanv2,
           aes(x = 1, y = Mean, fill = Group),
           width=1, color = "white") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = col_scheme2) +
  theme_void() +
  facet_nested_wrap(~Data + Method, ncol = 4,
                    nest_line = element_line(color = "grey60"),
                    solo_line = TRUE)+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0.5), unit = "cm"),
        panel.spacing = unit(1, "lines"),
        legend.position = "none")

pg <- ggplotGrob(P)
for(i in which(grepl("strip-t", pg$layout$name))){
  pg$grobs[[i]]$layout$clip <- "off"
}

grid::grid.draw(pg)

####cell type correlation####

peac <- read.csv("/media/gcpeac/Rachel/PEAC/Metadata/peacapril24.csv")
peac$SampleID_m0 <- trimws(peac$SampleID_m0)
peac$SampleID_m6 <- trimws(peac$SampleID_m6)

peacBL <- peac[peac$Visit == 1, ]

#generate meta

meta <- data.frame("Sample" = rownames(cellgeo$subclass$output))

load("/media/gcpeac/Rachel/PEAC/Metadata/PEAC_baseline_data.rdata")

meta$ID <- bld_metadata$QMULID[match(meta$Sample,
                                     bld_metadata$SAMID)]

meta$platelet <- peacBL$Platelets.V1[match(meta$ID,
                                           peacBL$SampleID_m0)]

meta$neutrophil <- peacBL$Neutrophils.V1[match(meta$ID,
                                               peacBL$SampleID_m0)]

meta$lymphocytes <- peacBL$Lymphocytes.V1[match(meta$ID,
                                                peacBL$SampleID_m0)]

add_data <- function(df){
  df$platelet <- meta$platelet[match(df$Patient,
                                     meta$Sample)]
  df$neutrophil <- meta$neutrophil[match(df$Patient,
                                         meta$Sample)]
  df$lymphocytes <- meta$lymphocytes[match(df$Patient,
                                           meta$Sample)]
  df
}

comb <- add_data(comb)

cell_plot <- function(df, subclass, cell, method, param = "Percentage"){
  sub <- df[df$Cluster == subclass & df$Method == method, ]
  sub$Cell <- sub[ , cell]
  sub$Param <- sub[ , param]
  
  stats <- cor.test(sub$Cell, sub$Param, method = "pearson")
  
  if(param == "Output"){
    ylab <- subclass
  }
  else{
    ylab <- paste(cell, "(%)")
  }
  
  sub <- sub[is.na(sub$Cell) == FALSE, ]
  
  ggplot(sub, aes(x = Cell, y = Param)) +
    geom_point(alpha = 0.8) +
    stat_smooth(method = "lm", col = "red") +
    labs(y = ylab, x = paste (cell, "(10^9/L)"),
         title = method,
         subtitle = paste0("r = ", formatC(stats$estimate,2),
                           " P = ", formatC(stats$p.value,2))) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 12),
          #plot.margin = margin(c(0, 0.5, 0, 0), unit = "cm")
          plot.margin = margin(c(0.1, 0.5, 0.1, 0.1), unit = "cm"))
}

music_tabula <- readRDS("music2_PEAC_tabula.rds")
music_tabula_percent <- as.data.frame(music_tabula$Est.prop.weighted * 100)
music_tabula_percent$Patient <- rownames(music_tabula_percent)
music_tabula_percent_long <- as.data.frame(pivot_longer(music_tabula_percent,
                                                        cols = -Patient,
                                                        names_to = "Cluster",
                                                        values_to = "Percentage"))
music_tabula_percent_long$Method <- "MuSiC"
music_tabula_percent_long <- add_data(music_tabula_percent_long)
music_tabula_percent_long$Cluster[music_tabula_percent_long$Cluster == "neutrophil"] <- "Neutrophil"

pdf("neutrophil_methods.pdf", width = 4.2, height = 2)
ggarrange(cell_plot(comb, "Neutrophil", "neutrophil", "CellGeometry") +
            scale_x_continuous(breaks = c(5, 10)),
          cell_plot(music_tabula_percent_long, "Neutrophil", "neutrophil", "MuSiC") +
            scale_x_continuous(breaks = c(5, 10)))
dev.off()

lym_plot <- function(df, groups = c("ILC", "B cells", "T cells"), method, param = "Percentage"){
  sub <- df[df$Group %in% groups & df$Method == method, ]
  sub$Param <- sub[ , param]
  
  sum <- sub %>% 
    group_by(Patient) %>%
    summarise("Sum" = sum(Param))
  
  sum$lymphocytes <- sub$lymphocytes[match(sum$Patient,
                                           sub$Patient)]
  
  stats <- cor.test(sum$lymphocytes, sum$Sum, method = "pearson")
  
  if(param == "Output"){
    ylab <- "Lymphocytes"
  }
  else{
    ylab <- "Lymphocytes (%)"
  }
  
  ggplot(sum, aes(x = lymphocytes, y = Sum)) +
    geom_point(alpha = 0.8) +
    stat_smooth(method = "lm", col = "red") +
    labs(y = ylab, x = "lymphocytes (10^9/L)",
         subtitle = paste0("r = ", formatC(stats$estimate,2),
                           " P = ", formatC(stats$p.value,2)),
         title = method) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 12),
          #plot.margin = margin(c(0, 0.5, 0, 0), unit = "cm")
          plot.margin = margin(c(0.1, 0.5, 0.1, 0.1), unit = "cm"))
}

cellgeo_plot <- list()
cellgeo_plot[["1"]] <- lym_plot(cellgeo_nest_percent_long, method = "CellGeometry") +
  scale_x_continuous(breaks = c(1, 2, 3))

plotlist <- lapply(c("MuSiC", "DWLS", "LinDeconSeq"), function(x){
  lym_plot(comb, method = x)+
    scale_x_continuous(breaks = c(1, 2, 3))
})

plotlist <- c(cellgeo_plot,
              plotlist)

pdf("PEAC_lymphocytes_method.pdf", width = 8.2, height = 2)
ggarrange(plotlist = plotlist, ncol = 4, nrow = 1)
dev.off()

####tabula####

temp <- readRDS("tabula_bld_markers_dualmeans.rds")

tabula_order <- c("B",
                  "Plasma",
                  "Plasmacytoid dendritic",
                  "Classical monocyte",
                  "Non-classical monocyte",
                  "Intermediate monocyte",
                  "Monocyte",
                  "Macrophage",
                  "Common myeloid progenitor",
                  "CD4-positive, alpha-beta T",
                  "CD8-positive, alpha-beta T",
                  "Naive thymus-derived CD4-positive, alpha-beta T",
                  "Regulatory T",
                  "Mature NK T",
                  "Natural killer",
                  "Basophil",
                  "Neutrophil",
                  "Erythrocyte",
                  "Platelet",
                  "Hematopoietic precursor",
                  "Hematopoietic stem")

tabula_col <- c("B" = "skyblue1",
                "Plasma" = "royalblue4",
                "Plasmacytoid dendritic" = "sienna1",
                "Classical monocyte" = "firebrick3",
                "Non-classical monocyte" = "firebrick1",
                "Intermediate monocyte" = "tomato1",
                "Monocyte" = "salmon1",
                "Macrophage" = "rosybrown1",
                "Common myeloid progenitor" = "burlywood", 
                "CD4-positive, alpha-beta T" = "purple3",
                "CD8-positive, alpha-beta T" = "plum4",
                "Naive thymus-derived CD4-positive, alpha-beta T",
                "Regulatory T" = "plum3",
                "Mature NK T" = "slateblue1",
                "Natural killer" = "aquamarine4",
                "Basophil" = "darkgoldenrod3",
                "Neutrophil" = "gold",
                "Erythrocyte" = "indianred",
                "Platelet" = "indianred4",
                "Hematopoietic precursor" = "brown1",
                "Hematopoietic stem" = "brown3")


bar_out <- function(df, percent, method, order_method){
  #Requires a nonnegative matrix.
  percent[is.na(percent) == TRUE] <- 0
  o <- seriate(percent, method = order_method)
  
  df <- df[df$Method == method, ]
  
  patient <- rownames(percent)[o[[1]]]
  
  df$Patient <- factor(df$Patient,
                       levels = patient)
  
  df$Percentage[is.na(df$Percentage) == TRUE] <- 0
  df <- df[df$Percentage > 0, ]
  
  df$Cluster <- as.vector(df$Cluster)
  df$Cluster <- gsub(" cell", "", df$Cluster)
  df$Cluster <- str_replace(df$Cluster, "^\\w{1}", toupper)
  df$Cluster <- factor(df$Cluster,
                       levels = tabula_order)
  
  
  ggplot(df) +
    geom_bar(aes(x = Patient, y = Percentage, fill = Cluster, color = Cluster), 
             position = "fill", stat = "identity") +
    scale_fill_manual(values = tabula_col, breaks = names(tabula_col), drop = FALSE)+
    scale_color_manual(values = tabula_col, breaks = names(tabula_col), drop = FALSE)+
    scale_y_continuous(labels = scales::percent_format(suffix = "")) +
    labs(x = "", y = "Cell (%)", title = method)+
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 8), #text size was 12 except for plot 14 and legend 11
          axis.text.x = element_blank(),
          axis.title = element_text(size = 8),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(fill=guide_legend(ncol=2)) 
}

music_tabula_percentv2 <- music_tabula_percent[ , -22]
NA_colsv3 <- colSums(music_tabula_percentv2) == 0
music_tabula_percentv2 <- music_tabula_percentv2[ , !NA_colsv3]

#stacked barplot
bar_out(music_tabula_percent_long, music_tabula_percentv2, "MuSiC", order_method)

####Table####

#detailing the mean cell proportions and its range 

library(dplyr)

stats_df <- comb %>% group_by(Cluster, Method) %>%
  summarise("Mean" = mean(Percentage),
            "Min" = min(Percentage),
            "Max" = max(Percentage))

stats_df$Mean <- format(round(stats_df$Mean, 2), nsmall = 2)
stats_df$Min <- format(round(stats_df$Min, 2), nsmall = 2)
stats_df$Max <- format(round(stats_df$Max, 2), nsmall = 2)

stats_df$Range <- paste0(trimws(stats_df$Mean),
                         " (",
                         trimws(stats_df$Min),
                         " - ",
                         trimws(stats_df$Max),
                         ")")

stats_df$Range[stats_df$Method == "MuSiC" & stats_df$Cluster %in% c("Basophil", "Neutrophil")] <- NA

stats_df_wide <- pivot_wider(stats_df,
                             id_cols = "Cluster",
                             names_from = "Method",
                             values_from = "Range")

stats_df_wide$Cluster <- gsub(" cells", "", stats_df_wide$Cluster)

# write.table(stats_df_wide, file = "prop_df.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

tabula_df <- music_tabula_percent_long %>% group_by(Cluster) %>%
  summarise("Mean" = mean(Percentage),
            "Min" = min(Percentage),
            "Max" = max(Percentage))

tabula_df$Mean <- format(round(tabula_df$Mean, 2), nsmall = 2)
tabula_df$Min <- format(round(tabula_df$Min, 2), nsmall = 2)
tabula_df$Max <- format(round(tabula_df$Max, 2), nsmall = 2)

tabula_df$Range <- paste0(trimws(tabula_df$Mean),
                          " (",
                          trimws(tabula_df$Min),
                          " - ",
                          trimws(tabula_df$Max),
                          ")")

tabula_df$Cluster <- gsub(" cell", "", tabula_df$Cluster)
tabula_df$Cluster <- str_replace(tabula_df$Cluster, "^\\w{1}", toupper)
tabula_df <- tabula_df[match(tabula_order, tabula_df$Cluster), ]

# write.table(tabula_df, file = "tabula_music_prop.txt",
#             append = FALSE, sep = "\t", col.names = NA, row.names = TRUE)


