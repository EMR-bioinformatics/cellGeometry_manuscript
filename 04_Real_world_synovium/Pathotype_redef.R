#Generating pathotype redefinition figures 

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/AMP_pathotype/Combine/cellgeo_outputs.rdata")

####B-M plot####
outputdf$Pathotype_orig <- as.vector(outputdf$Pathotype_orig)
outputdf$Pathotype_orig[is.na(outputdf$Pathotype_orig) == TRUE] <- "Ungraded"
outputdf$Pathotype_orig <- factor(outputdf$Pathotype_orig,
                                  levels = c("Fibroid",
                                             "Myeloid",
                                             "Lymphoid",
                                             "Ungraded"))

patho_col <- c("Lymphoid" = "steelblue",
               "Myeloid" = "firebrick",
               "Fibroid" = "green3",
               "Ungraded" = "grey80")

library(ggplot2)
library(ggpubr)
library(ggExtra)


BM_plot <- function(param, cohort, legend = FALSE){
  outputdf$Param <- outputdf[ , param]
  sub <- outputdf[outputdf$Cohort == cohort, ]
  
  if(grepl("Pathotype", param) == TRUE){
    P <- ggplot(sub, aes(x = B_cell, y = Monocyte, colour = Param)) +
      geom_point(alpha = 0.8, size = 1.3) +
      scale_color_manual(values = patho_col, name = "") +
      scale_y_continuous(trans='log10') +
      scale_x_continuous(trans='log10') +
      xlab("B cell") +
      ylab("Macrophage") +
      #ggtitle(cohort)+
      geom_vline(xintercept = 30, linetype = "dotted") +
      geom_hline(yintercept = 300, linetype = "dotted") +
      theme_classic() +
      theme(axis.text = element_text(color = "black", size = 8), #was 11
            axis.title = element_text(size = 10), #was 11
            plot.title = element_text(size = 12), #was 14
            legend.text = element_text(size = 11))
  } else{
    P <- ggplot(sub, aes(x = B_cell, y = Monocyte, colour = Param)) +
      geom_point(alpha = 0.8, size = 1.3) +
      scale_y_continuous(trans='log10') +
      scale_x_continuous(trans='log10') +
      scale_colour_gradientn(colours = paletteer_c("ggthemes::Classic Blue", 10),
                             na.value = "grey80",) +
      xlab("B cell") +
      ylab("Macrophage") +
      #ggtitle(cohort)+
      geom_vline(xintercept = 30, linetype = "dotted") +
      geom_hline(yintercept = 300, linetype = "dotted") +
      theme_classic() +
      theme(axis.text = element_text(color = "black", size = 8), #was 12
            axis.title = element_text( size = 10), #was 11
            plot.title = element_text( size = 12), #was 14
            legend.text = element_text(size = 12)) +
      labs(color = param)
  }
  
  if(legend == TRUE){
    P <- P + theme(legend.position = "top") #was "right"
  }
  else{
    P <- P + theme(legend.position = "none")
  }
  
  ggExtra::ggMarginal(P, type = "histogram")
}

param <- "Pathotype_orig"

ggarrange(BM_plot(param, "PEAC"),
          BM_plot(param, "STRAP"),
          BM_plot(param, "R4RA"),
          ncol = 1, nrow = 3)


####contingency tables####

####celltype status####

tab_df <- data.frame(table(outputdf$B_cell_status, outputdf$Monocyte_status, outputdf$Cohort))
colnames(tab_df)[colnames(tab_df) == "Var1"] <- "B_cell_status"
colnames(tab_df)[colnames(tab_df) == "Var2"] <- "Monocyte_status"
colnames(tab_df)[colnames(tab_df) == "Var3"] <- "Cohort"
tab_df$B_cell_status <- factor(tab_df$B_cell_status,
                               levels = c("Low", "High"))
tab_df$Monocyte_status <- factor(tab_df$Monocyte_status,
                                 levels = c("Low", "High"))
tab_df$Cohort <- factor(tab_df$Cohort,
                        levels = c("PEAC","STRAP", "R4RA"))

tab_df <- tab_df %>% group_by(Cohort) %>%
  mutate("Percentage" = Freq/sum(Freq) * 100)

tab_df$Percentage <- signif(tab_df$Percentage, 3)
tab_df$Annotate <- paste0(tab_df$Freq, "\n(", tab_df$Percentage, "%)")

ggplot(data = tab_df, aes(x = B_cell_status, y = Monocyte_status)) +
  geom_tile(fill = "white", color = "black") +
  geom_text(aes(label = Annotate)) +
  labs(x = "B cell status", y = "Macrophage status") +
  facet_wrap(~Cohort, ncol = 3) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(5, "lines"))

####original vs redefined status####

patho_df <- data.frame(table(outputdf$Pathotype_orig, outputdf$Pathotype_redef, outputdf$Cohort))
colnames(patho_df)[colnames(patho_df) == "Var1"] <- "Pathotype_orig"
colnames(patho_df)[colnames(patho_df) == "Var2"] <- "Pathotype_redef"
colnames(patho_df)[colnames(patho_df) == "Var3"] <- "Cohort"

patho_df <- patho_df %>% group_by(Cohort, Pathotype_redef) %>%
  mutate("Top" = ifelse(Freq == max(Freq), "Yes", "No"))
patho_df$Fill <- NA
patho_df$Fill[patho_df$Top == "Yes" & patho_df$Pathotype_redef == "Fibroid"] <- "#1B5E20FF"
patho_df$Fill[patho_df$Top == "No" & patho_df$Pathotype_redef == "Fibroid"] <- "#C8E6C9FF"

patho_df$Fill[patho_df$Top == "Yes" & patho_df$Pathotype_redef == "Myeloid"] <- "#B71C1CFF"
patho_df$Fill[patho_df$Top == "No" & patho_df$Pathotype_redef == "Myeloid"] <- "#FFCDD2FF"

patho_df$Fill[patho_df$Top == "Yes" & patho_df$Pathotype_redef == "Lymphoid"] <- "#0D47A1FF"
patho_df$Fill[patho_df$Top == "No" & patho_df$Pathotype_redef == "Lymphoid"] <- "#BBDEFBFF"

patho_df$Cohort <- factor(patho_df$Cohort,
                          levels = c("PEAC","STRAP", "R4RA"))

ggplot(data = patho_df, aes(x = Pathotype_redef, y = Pathotype_orig, fill = Fill)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq)) +
  scale_fill_identity() +
  labs(x = "Redefined pathotype", y = "Original pathotype") +
  facet_wrap(~Cohort, ncol = 3) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(5, "lines"))

####stacked pathotype####

library(tidyr)
library(dplyr)

arrange_data <- function(df, values, cohort){
  df$Sample <- rownames(df)
  test <- as.data.frame(pivot_longer(df,
                                     cols = -Sample,
                                     names_to = "Cluster",
                                     values_to = values))
  test$Cohort <- cohort
  
  test$Pathotype_orig <- outputdf$Pathotype_orig[match(test$Sample,
                                                       outputdf$Sample)]
  
  test$Pathotype_redef <- outputdf$Pathotype_redef[match(test$Sample,
                                                         outputdf$Sample)]
  
  test$Celltype <- NA
  test$Celltype[grepl("F", test$Cluster)] <- "Fibroblast"
  test$Celltype[grepl("M", test$Cluster)] <- "Macrophage"
  test$Celltype[grepl("B", test$Cluster)] <- "B cell"
  test$Celltype[grepl("T", test$Cluster)] <- "T cell"
  
  test$Celltype <- factor(test$Celltype,
                          levels = c("Fibroblast", "Macrophage",
                                     "T cell", "B cell"))
  test
}


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

#freq

peac_nest_output <- arrange_data(as.data.frame(peac1$nest_output), "Output", "PEAC")
r4ra_nest_output <- arrange_data(as.data.frame(r4ra$nest_output), "Output", "R4RA")
strap_nest_output <- arrange_data(as.data.frame(strap$nest_output), "Output", "STRAP")

comb_nest_output <- rbind(peac_nest_output,
                          r4ra_nest_output,
                          strap_nest_output)

patho_mean <- function(df, pathotype, param = "Percentage"){
  df$Param <- df[ , param]
  df$Patho <- df[ , pathotype]
  df <- df[df$Patho %in% c("Fibroid", "Myeloid", "Lymphoid"), ]
  df$Patho <- factor(df$Patho,
                     levels = c("Fibroid", "Myeloid", "Lymphoid"))
  mean <- as.data.frame(df %>% group_by(Cluster, Patho, Cohort) %>%
                          summarise("Mean" = mean(Param, na.rm = TRUE)))
  
  mean$Celltype <- df$Celltype[match(mean$Cluster,
                                     df$Cluster)]
  mean
  
}


patho_orig_mean <- patho_mean(comb_nest_output, "Pathotype_orig", "Output")
patho_redef_mean <- patho_mean(comb_nest_output, "Pathotype_redef", "Output")

patho_orig_mean$Definition <- "Original"
patho_redef_mean$Definition <- "Redefined"

patho_mean_comb <- rbind(patho_orig_mean,
                         patho_redef_mean)

patho_mean_comb$Cohort <- factor(patho_mean_comb$Cohort,
                                 levels = c("PEAC", "STRAP", "R4RA"))

library(ggplot2)
library(ggh4x)

ggplot() + 
  geom_col(data = patho_mean_comb, aes(x = Patho, y = Mean, fill = Cluster), 
           width = 0.5, color = "black") + 
  labs(x = "", y = "Cell freq") +
  scale_fill_manual(values = amp_scheme,
                    breaks = names(amp_scheme),
                    name = "",
                    guide = "none")+
  facet_nested(Celltype~Cohort + Definition, scales = "free_y", axes = "y",
               independent = "y", switch = "y")+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 13),
        plot.title = element_text(size = 14),
        strip.placement = "outside",
        panel.spacing.y = unit(1, "lines"),
        panel.spacing.x = unit(c(1,2,1,2,1), "lines"),
        strip.background = element_blank())


#percent

peac_nest_percent <- arrange_data(as.data.frame(peac1$nest_percent), "Percentage", "PEAC")
r4ra_nest_percent <- arrange_data(as.data.frame(r4ra$nest_percent), "Percentage", "R4RA")
strap_nest_percent <- arrange_data(as.data.frame(strap$nest_percent), "Percentage", "STRAP")

comb_nest_percent <- rbind(peac_nest_percent,
                           r4ra_nest_percent,
                           strap_nest_percent)

patho_orig_meanv2 <- patho_mean(comb_nest_percent, "Pathotype_orig")
patho_redef_meanv2 <- patho_mean(comb_nest_percent, "Pathotype_redef")

patho_orig_meanv2$Definition <- "Original"
patho_redef_meanv2$Definition <- "Redefined"

patho_mean_combv2 <- rbind(patho_orig_meanv2,
                           patho_redef_meanv2)

patho_mean_combv2$Cohort <- factor(patho_mean_combv2$Cohort,
                                   levels = c("PEAC", "STRAP", "R4RA"))


ggplot() + 
  geom_col(data = patho_mean_combv2, aes(x = Patho, y = Mean, fill = Cluster), 
           width = 0.5, color = "black") + 
  labs(x = "", y = "Cell (%)") +
  scale_fill_manual(values = amp_scheme,
                    breaks = names(amp_scheme),
                    name = "",
                    guide = "none")+
  facet_nested(Celltype~Cohort + Definition, scales = "free_y", axes = "y",
               independent = "y", switch = "y")+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 13),
        plot.title = element_text(size = 14),
        strip.placement = "outside",
        panel.spacing.y = unit(1, "lines"),
        panel.spacing.x = unit(c(1,2,1,2,1), "lines"),
        strip.background = element_blank())

####total cell counts####

peac_group_output <- arrange_data(as.data.frame(peac1$group$output), "Output", "PEAC")
r4ra_group_output <- arrange_data(as.data.frame(r4ra$group$output), "Output", "R4RA")
strap_group_output <- arrange_data(as.data.frame(strap$group$output), "Output", "STRAP")

comb_group_output <- rbind(peac_group_output,
                           r4ra_group_output,
                           strap_group_output)

patho_sum <- function(df, pathotype){
  df$Patho <- df[ , pathotype]
  df <- df[df$Patho %in% c("Fibroid", "Myeloid", "Lymphoid"), ]
  df$Patho <- factor(df$Patho,
                     levels = c("Fibroid", "Myeloid", "Lymphoid"))
  df$Sample <- factor(df$Sample)
  sum <- as.data.frame(df %>% group_by(Sample) %>%
                         summarise("Sum" = sum(Output, na.rm = TRUE)))
  sum$Patho <- df$Patho[match(sum$Sample,
                              df$Sample)]
  sum$Cohort <- df$Cohort[match(sum$Sample,
                                df$Sample)]
  as.data.frame(sum %>% group_by(Patho, Cohort) %>%
                  summarise("Mean" = mean(Sum, na.rm = TRUE),
                            "SD" = sd(Sum, na.rm = TRUE),
                            "num" = n(),
                            "SEM" = SD/sqrt(num)))
  
}

patho_orig_sum <- patho_sum(comb_group_output, "Pathotype_orig")
patho_redef_sum <- patho_sum(comb_group_output, "Pathotype_redef")

patho_orig_sum$Definition <- "Original"
patho_redef_sum$Definition <- "Redefined"

patho_sum_comb <- rbind(patho_orig_sum,
                        patho_redef_sum)

patho_sum_comb$Cohort <- factor(patho_sum_comb$Cohort,
                                levels = c("PEAC", "STRAP", "R4RA"))

ggplot() + 
  geom_errorbar(data = patho_sum_comb, aes(x = Patho, ymax = Mean + SEM,
                                           ymin = Mean - SEM), 
                width = 0.25)+
  geom_col(data = patho_sum_comb, aes(x = Patho, y = Mean, fill = Patho), 
           width = 0.5, color = "black") + 
  labs(x = "", y = "Cell freq") +
  scale_fill_manual(values = patho_col,
                    breaks = names(patho_col),
                    name = "",
                    guide = "none")+
  facet_nested(~Cohort + Definition, scales = "free", independent = "y") +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 13),
        plot.title = element_text(size = 14),
        panel.spacing.x = unit(c(1,2,1,2,1), "lines"),
        strip.background = element_blank())

####stacked barplot by pathotype####

sampledf <- outputdf[ , c("Sample", "Pathotype_redef", "Cohort")]
sampledf <- distinct(sampledf)

order_method = "CA"

cellgeo_order <- function(percent, patho, cohort){
  sub <- sampledf$Sample[sampledf$Pathotype_redef == patho & sampledf$Cohort == cohort]
  o <- seriate(percent[sub, ], method = order_method)
  patient <- rownames(percent[sub, ])
  patient[o[[1]]]
}

peac1_patient <- c(rev(cellgeo_order(as.data.frame(peac1$nest_percent), "Fibroid", "PEAC")),
                   rev(cellgeo_order(as.data.frame(peac1$nest_percent), "Myeloid", "PEAC")),
                   rev(cellgeo_order(as.data.frame(peac1$nest_percent), "Lymphoid", "PEAC")))

r4ra_patient <- c(rev(cellgeo_order(as.data.frame(r4ra$nest_percent), "Fibroid", "R4RA")),
                  rev(cellgeo_order(as.data.frame(r4ra$nest_percent), "Myeloid", "R4RA")),
                  rev(cellgeo_order(as.data.frame(r4ra$nest_percent), "Lymphoid", "R4RA")))

strap_patient <- c(cellgeo_order(as.data.frame(strap$nest_percent), "Fibroid", "STRAP"),
                   cellgeo_order(as.data.frame(strap$nest_percent), "Myeloid", "STRAP"),
                   cellgeo_order(as.data.frame(strap$nest_percent), "Lymphoid", "STRAP"))


patho_col <- c("Lymphoid" = "steelblue",
               "Myeloid" = "firebrick",
               "Fibroid" = "palegreen")
strip <- strip_themed(background_x = elem_list_rect(fill = alpha(rev(patho_col),
                                                                 0.6)))

comb_nest_percent$Pathotype_redef <- substring(comb_nest_percent$Pathotype_redef, 1, 1)
comb_nest_percent$Pathotype_redef <- factor(comb_nest_percent$Pathotype_redef, levels = c("F", "M", "L"))

bar_out <- function(patient, cohort){
  sub <- comb_nest_percent[comb_nest_percent$Cohort == cohort, ]
  sub$Sample <- factor(sub$Sample,
                       levels = patient)
  sub$Cluster <- factor(sub$Cluster,
                        levels = names(amp_scheme))
  
  ggplot(sub) +
    geom_col(aes(x = Sample, y = Percentage, fill = Cluster, color = Cluster), 
             position = "fill") +
    scale_fill_manual(values = amp_scheme, breaks = names(amp_scheme)) +
    scale_color_manual(values = amp_scheme, breaks = names(amp_scheme), guide = "none") +
    scale_y_continuous(labels = scales::percent_format(suffix = "")) +
    facet_grid2(~Pathotype_redef, scale = "free_x", space = "free", strip = strip)+
    coord_cartesian(clip = "off", ylim = c(0, 1)) +
    labs(x = "", y = "Cell (%)", title = cohort)+
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 8),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 8),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 10, vjust = 2),
          strip.background.x = element_rect(linewidth = 0),
          strip.text.x = element_text(margin = margin(b = 0)))
}

ggarrange(bar_out(peac1_patient, "PEAC"),
          bar_out(strap_patient, "STRAP"),
          bar_out(r4ra_patient, "R4RA"),
          ncol = 3, nrow = 1)






