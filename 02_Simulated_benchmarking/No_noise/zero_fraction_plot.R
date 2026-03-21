#visualising the outputs of the zero counts 

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/zero_fraction/simulated_dirichlet_zerofraction.rdata")

library(cellGeometry)
cellgeo <- readRDS("cellgeo_zerofraction.rds")
DWLS <- readRDS("DWLS_zerofraction.rds")
music <- readRDS("music_zerofraction.rds")
Lin <- readRDS("lin_zerofraction.rds")

# first retrieve all the results 

arrange_cellgeo <- function(list){
  
  df <- data.frame()
  for(i in names(list)){
    for(j in names(list[[i]])){
      for(x in names(list[[i]][[j]])){
        temp <- as.data.frame(list[[i]][[j]][[x]]$output$subclass$percent)
        
        temp$Sample <- rownames(temp)
        temp2 <- as.data.frame(pivot_longer(temp,
                                           cols = -Sample,
                                           names_to = "Cluster",
                                           values_to = "Pred"))
        temp2$zero_frac <- i
        temp2$Times <- as.numeric(gsub("Times ", "", j))
        temp2$Rep <- as.numeric(gsub("Rep ", "", x))
        temp2$Method <- "CellGeometry"
        df <- rbind(df,
                    temp2)
      }
    }
  }
  
  df
}

arrange_music <- function(list){
  
  df <- data.frame()
  for(i in names(list)){
    for(j in names(list[[i]])){
      for(x in names(list[[i]][[j]])){
        temp <- as.data.frame(list[[i]][[j]][[x]]$output$Est.prop.weighted) * 100
        
        temp$Sample <- rownames(temp)
        temp2 <- as.data.frame(pivot_longer(temp,
                                            cols = -Sample,
                                            names_to = "Cluster",
                                            values_to = "Pred"))
        temp2$zero_frac <- i
        temp2$Times <- as.numeric(gsub("Times ", "", j))
        temp2$Rep <- as.numeric(gsub("Rep ", "", x))
        temp2$Method <- "MuSiC"
        df <- rbind(df,
                    temp2)
      }
    }
  }
  
  df
}

arrange_DWLS <- function(list){
  
  df <- data.frame()
  for(i in names(list)){
    for(j in names(list[[i]])){
      for(x in names(list[[i]][[j]])){
        temp <- as.data.frame(t(list[[i]][[j]][[x]]$output)) * 100
        
        temp$Sample <- rownames(temp)
        temp2 <- as.data.frame(pivot_longer(temp,
                                            cols = -Sample,
                                            names_to = "Cluster",
                                            values_to = "Pred"))
        temp2$zero_frac <- i
        temp2$Times <- as.numeric(gsub("Times ", "", j))
        temp2$Rep <- as.numeric(gsub("Rep ", "", x))
        temp2$Method <- "DWLS"
        df <- rbind(df,
                    temp2)
      }
    }
  }
  
  df
}

arrange_lin <- function(list){
  
  df <- data.frame()
  for(i in names(list)){
    for(j in names(list[[i]])){
      for(x in names(list[[i]][[j]])){
        temp <- as.data.frame(list[[i]][[j]][[x]]$output) * 100
        
        temp$Sample <- rownames(temp)
        temp2 <- as.data.frame(pivot_longer(temp,
                                            cols = -Sample,
                                            names_to = "Cluster",
                                            values_to = "Pred"))
        temp2$zero_frac <- i
        temp2$Times <- as.numeric(gsub("Times ", "", j))
        temp2$Rep <- as.numeric(gsub("Rep ", "", x))
        temp2$Method <- "LinDeconSeq"
        df <- rbind(df,
                    temp2)
      }
    }
  }
  
  df
}



comb <- rbind(arrange_cellgeo(cellgeo),
              arrange_music(music),
              arrange_DWLS(DWLS),
              arrange_lin(Lin))

comb$Cluster <- gsub("\\.", "-", comb$Cluster)

#get actual proportions 

sim_percent_long <- as.data.frame(sim_percent)
sim_percent_long$Patient <- rownames(sim_percent_long)

sim_percent_long <-  as.data.frame(pivot_longer(sim_percent_long,
                                               cols = -Patient,
                                               names_to = "Cluster",
                                               values_to = "Percentage"))

sim_percent_long20 <- as.data.frame(sim_percent20)
sim_percent_long20$Patient <- rownames(sim_percent_long20)

sim_percent_long20 <-  as.data.frame(pivot_longer(sim_percent_long20,
                                                cols = -Patient,
                                                names_to = "Cluster",
                                                values_to = "Percentage"))


comb$Obs <- NA
comb$Obs[comb$zero_frac == "zero_10"] <- sim_percent_long$Percentage[match(paste(comb$Sample[comb$zero_frac == "zero_10"],
                                                                                 comb$Cluster[comb$zero_frac == "zero_10"]),
                                                                           paste(sim_percent_long$Patient,
                                                                                 sim_percent_long$Cluster))]

comb$Obs[comb$zero_frac == "zero_20"] <- sim_percent_long20$Percentage[match(paste(comb$Sample[comb$zero_frac == "zero_20"],
                                                                                 comb$Cluster[comb$zero_frac == "zero_20"]),
                                                                           paste(sim_percent_long20$Patient,
                                                                                 sim_percent_long20$Cluster))]
#add the other zero_frac

add_sim <- function(sim, zero_frac){
  df_long <- as.data.frame(sim)
  df_long$Patient <- rownames(df_long)
  
  df_long <-  as.data.frame(pivot_longer(df_long,
                                        cols = -Patient,
                                        names_to = "Cluster",
                                        values_to = "Percentage"))
  
  comb$Obs[comb$zero_frac == zero_frac] <- df_long$Percentage[match(paste(comb$Sample[comb$zero_frac == zero_frac],
                                                                          comb$Cluster[comb$zero_frac == zero_frac]),
                                                                    paste(df_long$Patient,
                                                                          df_long$Cluster))]
  
  comb
}


comb <- add_sim(sim_percent_none, "zero_0")
comb <- add_sim(sim_percent5, "zero_5")
comb <- add_sim(sim_percent50, "zero_50")

####obs/pred plot####

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

mk <- cellgeo$zero_10$`Times 30`$`Rep 1`$output$mk

# music_percent <- as.data.frame(music$zero_10$`Times 30`$`Rep 1`$output$Est.prop.weighted) * 100
# music_percent <- music_percent[ , match(colnames(sim_percent),
#                                         colnames(music_percent))]
# DWLS_percent <- as.data.frame(t(DWLS$zero_10$`Times 30`$`Rep 1`$output)) * 100
# DWLS_percent <- DWLS_percent[ , match(colnames(sim_percent),
#                                       colnames(DWLS_percent))]
# 
# Lin_percent <- as.data.frame(Lin$zero_10$`Times 30`$`Rep 1`$output) * 100
# colnames(Lin_percent) <- gsub("\\.", "-", colnames(Lin_percent))
# Lin_percent <- Lin_percent[ , match(colnames(sim_percent),
#                                     colnames(Lin_percent))]
# 
# plot_pred(sim_percent, cellgeo$zero_10$`Times 30`$`Rep 1`$output$subclass$percent,
#           scheme = amp_scheme)
# plot_pred(sim_percent, as.matrix(music_percent), mk, scheme = amp_scheme)
# plot_pred(sim_percent, as.matrix(DWLS_percent), mk, scheme = amp_scheme)
# plot_pred(sim_percent, as.matrix(Lin_percent), mk, scheme = amp_scheme)

# rm(DWLS_percent)
# rm(Lin_percent)
# rm(music_percent)

scatter_out <- function(sim, zero_frac){
  music_percent <- music[[zero_frac]]$`Times 30`$`Rep 1`$output$Est.prop.weighted * 100
  music_percent <- music_percent[ , match(colnames(sim),
                                          colnames(music_percent))]

  DWLS_percent <- t(DWLS[[zero_frac]]$`Times 30`$`Rep 1`$output) * 100
  DWLS_percent <- DWLS_percent[ , match(colnames(sim),
                                        colnames(DWLS_percent))]

  Lin_percent <- Lin[[zero_frac]]$`Times 30`$`Rep 1`$output * 100
  colnames(Lin_percent) <- gsub("\\.", "-", colnames(Lin_percent))
  Lin_percent <- Lin_percent[ , match(colnames(sim),
                                      colnames(Lin_percent))]
  
  cellgeo_rsq <- format(Rsq(sim, cellgeo[[zero_frac]]$`Times 30`$`Rep 1`$output$subclass$percent), digits = 3)
  music_rsq <- format(Rsq(sim, music_percent), digits = 3)
  dwls_rsq <- format(Rsq(sim, DWLS_percent), digits = 3)
  lin_rsq <- format(Rsq(sim, as.matrix(Lin_percent)), digits = 3)
  
  fig <- ggarrange(plot_pred(sim, cellgeo[[zero_frac]]$`Times 30`$`Rep 1`$output$subclass$percent,
             scheme = amp_scheme) + labs(title = "CellGeometry", subtitle =  bquote(R^2 ~"="~ .(cellgeo_rsq))) +
               theme(plot.title = element_text(margin=margin(0,0,0.2,0))),
            plot_pred(sim, music_percent, mk, scheme = amp_scheme) + labs(title = "MuSiC", subtitle =  bquote(R^2 ~"="~ .(music_rsq)))+
              theme(plot.title = element_text(margin=margin(0,0,0.2,0))),
            plot_pred(sim, DWLS_percent, mk, scheme = amp_scheme) + labs(title = "DWLS", subtitle =  bquote(R^2 ~"="~ .(dwls_rsq)))+
              theme(plot.title = element_text(margin=margin(0,0,0.2,0))),
            plot_pred(sim, as.matrix(Lin_percent), mk, scheme = amp_scheme) + labs(title = "LinDeconSeq", subtitle =  bquote(R^2 ~"="~ .(lin_rsq)))+
              theme(plot.title = element_text(margin=margin(0,0,0.2,0))),
            ncol = 4)
  
  annotate_figure(fig, top = text_grob(paste0(gsub(".*_", "", zero_frac), "% subclasses set to 0 in each sample"),
                                       size = 10, color = "black", face = "bold"))
}

Rsq <- function(obs, pred) {
  if (anyNA(pred)) {
    pred[is.na(pred)] <- 0
  }
  rss <- sum((pred - obs)^2)
  tss <- sum((obs - mean(obs))^2)
  1 - rss/tss
}

library(ggpubr)
library(ggplot2)

ggarrange(scatter_out(sim_percent_none, "zero_0"),
          scatter_out(sim_percent5, "zero_5"),
          scatter_out(sim_percent, "zero_10"),
          scatter_out(sim_percent20, "zero_20"),
          scatter_out(sim_percent50, "zero_50"),
          nrow = 5)

####failure rate of DWLS####

failure <- data.frame()

for(i in names(DWLS)){
  for(j in names(DWLS[[i]])){
    for(x in names(DWLS[[i]][[j]])){
    samples <- sum(colSums(is.na(DWLS[[i]][[j]][[x]]$output)) > 0)
    temp <- data.frame("Zero_frac" = i,
                       "Rep" = as.numeric(gsub("Rep ", "", x)),
                       "Times" =  as.numeric(gsub("Times ", "", j)),
                       "Samples" = samples,
                       "Percentage" = (samples/25)*100)
    failure <- rbind(failure,
                     temp)
    }
  }
}

failure_range <- failure %>%
  group_by(Zero_frac) %>%
  summarise("Min" = min(Percentage),
            "Max" = max(Percentage),
            "Mean" = mean(Percentage))

failure_range$Zero_frac <- as.numeric(gsub(".*_", "", failure_range$Zero_frac))

ggplot(failure_range, aes(x = Mean, y = Zero_frac)) +
  geom_point() +
  geom_errorbar(aes(xmax = Max, xmin = Min), width = 1.2) +
  labs(title = "DWLS",
       x = "Quadratic programming failure rate (%)",
       y = "% subclasses set to 0") +
  xlim(c(0, max(failure_range$Max))) +
  scale_y_continuous(breaks = c(0, 5, 10, 20, 50)) +
  coord_flip()+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10), #was 8
        axis.title = element_text(size = 11), # was 10
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size = 10))

####Rsq mean + SEM with zero_frac####

rsq_df <- data.frame()

for(i in levels(comb$Method)){
  for(x in unique(comb$zero_frac)){
    for(j in unique(comb$Rep)){
    temp <- data.frame("R2" = Rsq(comb$Obs[comb$Method == i & comb$zero_frac == x & comb$Rep == j], 
                                  comb$Pred[comb$Method == i & comb$zero_frac == x & comb$Rep == j]),
                       "Method" = i,
                       "zero_frac" = x,
                       "Rep" = j)
    rsq_df <- rbind(rsq_df,
                    temp)
    }
  }
}

rsq_df$zero_frac <- as.numeric(gsub("zero_", "", rsq_df$zero_frac))
rsq_df$Method <- factor(rsq_df$Method,
                        levels = c("CellGeometry",
                                   "MuSiC",
                                   "DWLS",
                                   "LinDeconSeq"))

rsq_mean <- rsq_df %>%
  group_by(Method, zero_frac) %>%
  summarise(Mean = mean(R2),
            SD = sd(R2),
            N = n(),
            SEM = SD/sqrt(N))


ggplot() +
  geom_smooth(data = rsq_df, aes(x = zero_frac, y = R2,
                                 color = Method, group = Method),
              method = "loess", alpha = 0.2, linewidth = 0.5, se = FALSE,
              show.legend = TRUE) +
  geom_errorbar(data = rsq_mean,
                aes(x = zero_frac, y = Mean, color = Method, group = Method,
                    ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.8,
                linewidth = 0.4, 
                show.legend = FALSE)  +
  geom_point(data = rsq_mean,
             aes(x = zero_frac, y = Mean, fill = Method, 
                 group = Method), pch = 21,
             show.legend = TRUE) + 
  labs(x = "% of subclasses set to 0", y = "Rsq") +
  scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "LinDeconSeq" = "#00B0F6",
                                "DWLS" = "#E76BF3"),
                     name = "method",
                     drop = FALSE) +
  scale_fill_manual(values = c("CellGeometry" = "#F8766D",
                                "MuSiC" = "#00BF7D",
                                "LinDeconSeq" = "#00B0F6",
                                "DWLS" = "#E76BF3"),
                     name = "method",
                     drop = FALSE) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10), 
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 11), 
        plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"))



