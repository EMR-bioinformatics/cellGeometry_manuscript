####Background####

# plotting the metrics of the noise simulated deconvolution

setwd("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/Noise_updated")

sim_sampled_dir_all <- readRDS("../simulated_dirichlet_all.rds")
load("../simulated_dirichlet.rdata")

####Metric plot####

comb_all_noise <- readRDS("music_lin_noise_out.rds")
temp <- readRDS("cellgeo_DWLS_noise_out.rds")

comb_all_noise <- rbind(comb_all_noise,
                        temp)

rm(temp)

#remove previous sqrt noise and replace with updated sqrt noise

comb_all_noise <- comb_all_noise[comb_all_noise$Type != "Sqrt noise", ]

cellgeo_sqrt <- readRDS("cellgeo_sqrt_noise_correct.rds")
DWLS_sqrt <- readRDS("DWLS_sqrt_noise_correct.rds")
music_sqrt <- readRDS("music_sqrt_noise_correct.rds")
lin_sqrt <- readRDS("lin_sqrt_noise_correct.rds")

percent_retrieve <- function(list){
  df <- data.frame()
  for(i in names(list)){
    for(x in names(list[[i]])){
      temp <- as.data.frame(list[[i]][[x]]$metrics_percent)
      temp$Rep <- as.numeric(gsub("Rep ", "", x))
      temp$Noise <- i
      temp <- add_column(temp,
                         "Cluster" = rownames(temp),
                         .before = 1)
      df <- rbind(df,
                  temp)
    }
  }
  df
}

library(tibble)

cellgeo_baseline <- readRDS("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Benchmarking/AMP/Dirichlet/cellgeo_dirichlet_output500.rds")
DWLS_baseline <- readRDS("../DWLS_dirichlet_output.rds")
music_baseline <- readRDS("../music_dirichlet_output.rds")
lin_baseline <- readRDS("../Lin_dirichlet_output.rds")

cellgeo_sqrt_noise <- rbind(percent_retrieve(cellgeo_baseline),
                            percent_retrieve(cellgeo_sqrt))
cellgeo_sqrt_noise$Noise[cellgeo_sqrt_noise$Noise == "Times 30"] <- "SD = 0"
cellgeo_sqrt_noise$Method <- "CellGeometry"

DWLS_sqrt$`SD = 0` <- DWLS_baseline$`Times 30`
DWLS_sqrt_noise <- percent_retrieve(DWLS_sqrt)
DWLS_sqrt_noise$Method <- "DWLS"

music_sqrt$`SD = 0` <- music_baseline$`Times 30`
music_sqrt_noise <- percent_retrieve(music_sqrt)
music_sqrt_noise$Method <- "MuSiC"

lin_sqrt$`SD = 0` <- lin_baseline$`Times 30`
lin_sqrt_noise <- percent_retrieve(lin_sqrt)
lin_sqrt_noise$Method <- "LinDeconSeq"

sqrt_noise <- rbind(cellgeo_sqrt_noise,
                    DWLS_sqrt_noise,
                    music_sqrt_noise,
                    lin_sqrt_noise)

sqrt_noise$Type <- "Sqrt noise"
sqrt_noise$Noise <- gsub("SD = ", "", sqrt_noise$Noise)

comb_all_noise <- rbind(comb_all_noise,
                        sqrt_noise)

npass <- readRDS("cellgeo_shift_npass_var2point5.rds")

npass_df <- percent_retrieve(npass)
npass_df$Method <- "CellGeometry 2-npass"
npass_df$Type <- "Shift noise"

npass_df$Noise <- gsub("SD = ", "", npass_df$Noise)

comb_all_noise <- rbind(comb_all_noise,
                        npass_df)

comb_all_noise$Method <- factor(comb_all_noise$Method,
                                levels = c("CellGeometry",
                                           "CellGeometry 2-npass",
                                           "MuSiC",
                                           "DWLS",
                                           "LinDeconSeq"))


comb_meandf <- comb_all_noise %>% 
  group_by(Noise, Method, Type) %>%
  summarise(Mean = mean(RMSE),
            SD = sd(RMSE),
            N = n(),
            SEM = Mean/sqrt(N))

library(ggplot2)
library(ggnewscale)
library(dplyr)

metric_plot <- function(df, metric, lim){
  df$Metric <- df[ , metric]
  
  meandf <- df %>% 
    group_by(Noise, Method) %>%
    summarise(Mean = mean(Metric),
              SD = sd(Metric),
              N = n(),
              SEM = Mean/sqrt(N))
  
  ggplot() +
    geom_smooth(data = df, se = FALSE,
                aes(x = as.numeric(Noise), y = Metric,
                    color = Method, group = Method),
                method = "loess", alpha = 0.2, linewidth = 0.5, 
                show.legend = TRUE) +
    scale_color_manual(values = c("CellGeometry" = "#F8766D",
                                  "CellGeometry 2-npass" = "#c25048",
                                  "MuSiC" = "#00BF7D",
                                  "LinDeconSeq" = "#00B0F6",
                                  "DWLS" = "#E76BF3"),
                       name = "method",
                       drop = FALSE) +
    new_scale_color() +
    geom_errorbar(data = meandf,
                  aes(x = as.numeric(Noise), y = Mean, color = Method, group = Method,
                      ymin = Mean - SEM, ymax = Mean + SEM),
                  width = diff(range(as.numeric(meandf$Noise))) *0.02,
                  linewidth = 0.4, 
                  show.legend = FALSE)  +
    geom_point(data = meandf,
               aes(x = as.numeric(Noise), y = Mean, fill = Method, color = Method,
                   group = Method,
                   shape = Method),
               show.legend = TRUE) + #, size = 1.6) +
    labs(x = "noise", y = metric) +
    scale_color_manual(values = c("CellGeometry" = "black",
                                    "CellGeometry 2-npass" = "#c25048",
                                    "MuSiC" = "#00BF7D",
                                    "LinDeconSeq" = "#00B0F6",
                                    "DWLS" = "#E76BF3"),
                         guide = FALSE,
                         drop = FALSE) +
    scale_shape_manual(values = c("CellGeometry" = 24,
                                  "CellGeometry 2-npass" = 21,
                                  "MuSiC" = 21,
                                  "LinDeconSeq" = 21,
                                  "DWLS" = 21),
                       guide = FALSE)+
    scale_fill_manual(values = c("CellGeometry" = "#F8766D",
                                  "CellGeometry 2-npass" = "#c25048",
                                  "MuSiC" = "#00BF7D",
                                  "LinDeconSeq" = "#00B0F6",
                                  "DWLS" = "#E76BF3"),
                       name = "method",
                       drop = FALSE) +
    coord_cartesian(ylim = lim) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 9), # was 8
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(size = 10), #was 9
          plot.title = element_text(size = 9),
          legend.text = element_text(size = 8), #was 9
          legend.title = element_text(size = 10), #was 9
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm")) +
          #plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    guides(fill =  guide_legend(override.aes = list(shape = c(24, rep(21, 4)),
                                                    color = c("CellGeometry" = "black",
                                                              "CellGeometry 2-npass" = "#c25048",
                                                              "MuSiC" = "#00BF7D",
                                                              "DWLS" = "#E76BF3",
                                                              "LinDeconSeq" = "#00B0F6"))))
}

library(ggpubr)

pdf("RMSE_allv2_npassv2_correct.pdf", width = 12.5, height = 3, onefile=FALSE)
ggarrange(metric_plot(comb_all_noise[comb_all_noise$Type == "Gaussian noise" &
                                       comb_all_noise$Noise %in% c(0, 100, 200,
                                                                   300, 400, 500), ],
                      "RMSE", c(0.5, 7)) + ggtitle("Gaussian noise") ,
          metric_plot(comb_all_noise[comb_all_noise$Type == "Log noise", ],
                      "RMSE", c(0.5, 7)) + ggtitle("log noise"),
          metric_plot(comb_all_noise[comb_all_noise$Type == "Sqrt noise" &
                                       comb_all_noise$Noise %in% c(0, 10, 25, 
                                                                   50, 75, 100), ],
                      "RMSE", c(0.5, 7)) + ggtitle("sqrt noise"),
          metric_plot(comb_all_noise[comb_all_noise$Type == "Shift noise", ],
                      "RMSE", c(0.5, 7)) + ggtitle("shift noise"),
          ncol = 4, common.legend = TRUE, legend = "right")
dev.off()

pdf("Rsq_allv2_npassv2_correct.pdf", width = 12.5, height = 3, onefile=FALSE)
ggarrange(metric_plot(comb_all_noise[comb_all_noise$Type == "Gaussian noise" &
                                       comb_all_noise$Noise %in% c(0, 100, 200,
                                                                   300, 400, 500), ],
                      "Rsq", c(-8, 1)) + ggtitle("Gaussian noise"),
          metric_plot(comb_all_noise[comb_all_noise$Type == "Log noise", ],
                      "Rsq", c(-8, 1)) + ggtitle("Log noise"),
          metric_plot(comb_all_noise[comb_all_noise$Type == "Sqrt noise" &
                                       comb_all_noise$Noise %in% c(0, 10, 25, 
                                                                   50, 75, 100), ],
                      "Rsq", c(-8, 1)) + ggtitle("Sqrt noise"),
          metric_plot(comb_all_noise[comb_all_noise$Type == "Shift noise", ],
                      "Rsq", c(-8, 1)) + ggtitle("Shift noise"),
          ncol = 4, common.legend = TRUE, legend = "right")
dev.off()

####noise plots####

noise_all <- readRDS("default_noise.rds")
shift_all <- readRDS("shift_noise.rds")
add_all <- readRDS("add_noise.rds")
log_all <- readRDS("log_noise.rds")
sqrt_all <- readRDS("sqrt_noise_correct.rds")

library(ggrastr)

noise_plot <- function(orig, noise, title, alpha = 0.05){
  plotdf <- data.frame("Original" = orig[ , 1],
                       "Noise" = noise[ , 1])
  
  ggplot(plotdf, aes(x = Original, y = Noise)) +
    rasterize(geom_point(alpha = alpha, size = 0.75), scale = 0.6,
               #dpi = 2000 #for main) +
              dpi = 500)+
    theme_classic()+
    scale_x_log10()+
    scale_y_log10()+
    ggtitle(title) + 
    theme(axis.text = element_text(color = "black", size = 8),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 10))
  
}

orig <- sim_sampled_dir_all$`Times 30`$`Rep 1`

pdf("example_noise_SD_log_correct.pdf", 
    width = 10, height = 2.5)
ggarrange(noise_plot(orig, noise_all$add$`Rep 1`, "Gaussian noise"),
          noise_plot(orig, log_all$`SD = 0.25`$`Rep 1`, "Log noise"),
          noise_plot(orig, sqrt_all$`SD = 10`$`Rep 1`, "Sqrt noise"),
          noise_plot(orig, noise_all$shift$`Rep 1`, "Shift noise"),
          ncol = 4)
dev.off()

pdf("add_noise_SD_logv2_correct.pdf", 
width = 12.5, height = 2.5)
ggarrange(noise_plot(orig, noise_all$add$`Rep 1`, "SD = 100"),
          noise_plot(orig, add_all$`SD = 250`$`Rep 1`, "SD = 200"),
          noise_plot(orig, add_all$`SD = 500`$`Rep 1`, "SD = 300"),
          noise_plot(orig, add_all$`SD = 750`$`Rep 1`, "SD = 400"),
          noise_plot(orig, add_all$`SD = 1000`$`Rep 1`, "SD = 500"),
          ncol = 5)
dev.off()

pdf("log_noise_SD_log_correct.pdf", 
    width = 12.5, height = 2.5)
ggarrange(noise_plot(orig, noise_all$log$`Rep 1`, "SD = 0.1"),
          noise_plot(orig, log_all$`SD = 0.25`$`Rep 1`, "SD = 0.25"),
          noise_plot(orig, log_all$`SD = 0.5`$`Rep 1`, "SD = 0.5"),
          noise_plot(orig, log_all$`SD = 0.75`$`Rep 1`, "SD = 0.75"),
          noise_plot(orig, log_all$`SD = 1`$`Rep 1`, "SD = 1"),
          ncol = 5)
dev.off()

pdf("sqrt_noise_SD_logv2_correct.pdf", 
    width = 12.5, height = 2.5)
ggarrange(noise_plot(orig, sqrt_all$`SD = 10`$`Rep 1`, "SD = 10"),
          noise_plot(orig, sqrt_all$`SD = 25`$`Rep 1`, "SD = 25"),
          noise_plot(orig, sqrt_all$`SD = 50`$`Rep 1`, "SD = 50"),
          noise_plot(orig, sqrt_all$`SD = 75`$`Rep 1`, "SD = 75"),
          noise_plot(orig, sqrt_all$`SD = 100`$`Rep 1`, "SD = 100"),
          ncol = 5)
dev.off()

pdf("shift_noise_SD_log_correct.pdf", 
    width = 12.5, height = 2.5)
ggarrange(noise_plot(orig, noise_all$shift$`Rep 1`, "SD = 0.5"),
          noise_plot(orig, shift_all$`SD = 0.75`$`Rep 1`, "SD = 0.75"),
          noise_plot(orig, shift_all$`SD = 1`$`Rep 1`, "SD = 1"),
          noise_plot(orig, shift_all$`SD = 1.5`$`Rep 1`, "SD = 1.5"),
          noise_plot(orig, shift_all$`SD = 2`$`Rep 1`, "SD = 2"),
          ncol = 5)
dev.off()

####DWLS sample completion####

completion_out <- function(DWLS, type){
  
  completion <- data.frame()
  
  for(i in names(DWLS)){
    for(x in names(DWLS[[i]])){
      samples <- sum(colSums(is.na(DWLS[[i]][[x]]$output)) > 0)
      temp <- data.frame("Noise" = gsub("SD = ", "", i),
                         "Rep" = gsub("Rep ", "", x),
                         "Method" = "DWLS",
                         "Samples" = 25-samples,
                         "Percentage" = ((25 -samples)/25)*100)
      completion <- rbind(completion,
                          temp)
    }
  }
  
  completion$Rep <- factor(completion$Rep)
  completion$Noise <- factor(as.numeric(completion$Noise))
  completion$Type <- type
  completion
}

DWLS_add <- readRDS("DWLS_add_noise.rds")
DWLS_add$`SD = 0` <- DWLS_baseline$`Times 30`
DWLS_log <- readRDS("DWLS_log_noise.rds")
DWLS_log$`SD = 0` <- DWLS_baseline$`Times 30`
DWLS_shift <- readRDS("DWLS_shift_noise.rds")
DWLS_shift$`SD = 0` <- DWLS_baseline$`Times 30`

add_completion <- completion_out(DWLS_add, "Gaussian noise")
log_completion <- completion_out(DWLS_log, "Log noise")
sqrt_completion <- completion_out(DWLS_sqrt, "Sqrt noise")
shift_completion <- completion_out(DWLS_shift, "Shift noise")

comb_completion <- rbind(add_completion[add_completion$Noise %in% c(0, 100, 200,
                                                                      300, 400, 500), ],
                         log_completion,
                         sqrt_completion[sqrt_completion$Noise %in% c(0, 10, 25,
                                                                        50, 75, 100), ],
                         shift_completion)

comb_completion$Noise <- as.numeric(as.vector(comb_completion$Noise))
comb_completion$Type <- factor(comb_completion$Type,
                               levels = c("Gaussian noise",
                                          "Log noise",
                                          "Sqrt noise",
                                          "Shift noise"))

comb_completion$Noise <- factor(comb_completion$Noise,
                                levels = sort(unique(comb_completion$Noise)))

comb_completion_range <- comb_completion %>%
  group_by(Type, Noise) %>%
  summarise("Min" = min(Percentage),
            "Max" = max(Percentage),
            "Mean" = mean(Percentage))

comb_completion_range$Noise <- factor(comb_completion_range$Noise,
                                      levels = rev(sort(unique(comb_completion$Noise))))

pdf("DWLS_log_shift_failure.pdf", width = 5.5, height = 3)
ggplot(comb_completion_range[comb_completion_range$Type %in% c("Log noise", "Shift noise"), ]) +
  geom_point(aes(x = Noise, y = 100 - Mean)) +
  geom_errorbar(aes(x = Noise, ymax = 100 - Max, ymin = 100 - Min),
                width = 0.2) +
  theme_classic() +
  labs(y = "Quadratic programming failure rate (%)",
       x = "noise") +
  ylim(c(0, 100)) +
  coord_flip() +
  facet_wrap(~Type, scales = "free_y") +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 10, hjust = 0.5))
dev.off()

