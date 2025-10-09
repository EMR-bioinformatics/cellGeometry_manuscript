
# human brain atlas 1.0
# https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443

# all neurons
# downloaded
# https://datasets.cellxgene.cziscience.com/c2f66cd5-4ff4-4578-876c-55783a57cf8f.h5ad
# 33 Gb
# 2,480,956 cells

# all non-neuron cells
# 888,263 cells
# https://datasets.cellxgene.cziscience.com/99f27be8-9fac-451e-9723-9e4c7191589e.h5ad
# 4 Gb

library(zellkonverter)
library(SingleCellExperiment)

library(devtools)
setwd("/users/myles/documents/github/cellGeometry")
load_all()

library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

brain <- readH5AD("/Users/myles/R/Deconv/c2f66cd5-4ff4-4578-876c-55783a57cf8f.h5ad",
                  use_hdf5 = TRUE, reader = "R")

mat <- brain@assays@data$X
rownames(mat) <- rownames(brain)  # need to add rownames (genes)
meta <- brain@colData@listData

str(mat)
sort(table(meta$roi))
sort(table(meta$supercluster_term))

# other possible columns:
# dissection - fine subregions
# ROIgroupfine - similar to tissue region

mk <- cellMarkers(mat, subclass = meta$roi,
                  cellgroup = meta$supercluster_term, 
                  dual_mean = TRUE, cores = 8)
# 12.2 + 14.7 mins (intel 8 cores)

mk <- gene2symbol(mk, ensDb_v110)

signature_heatmap(mk, show_row_names = F, row_title_rot = 0, column_title_rot = 45)

# non-neuronal cells
brainNN <- readH5AD("/Users/myles/R/Deconv/99f27be8-9fac-451e-9723-9e4c7191589e.h5ad",
                  use_hdf5 = TRUE, reader = "R")

mat2 <- brainNN@assays@data$X
rownames(mat2) <- rownames(brainNN)  # need to add rownames (genes)
meta2 <- brainNN@colData@listData

sort(table(meta2$supercluster_term))
sort(table(meta2$cell_type))

mkNN <- cellMarkers(mat2, subclass = meta2$cell_type,
                    cellgroup = meta2$supercluster_term,
                    dual_mean = TRUE, cores = 8)
# 4.29 + 4.18 mins (intel 8 cores)

mkNN <- gene2symbol(mkNN, ensDb_v110)

# mkm <- mergeMarkers(mk, mkNN)
# plot(mkm$qqmerge)
mkm <- mergeMarkers(mk, mkNN, transform = "none")
mkm <- updateMarkers(mkm, expfilter = 0.2, nsubclass = 25)

# save(file = "/Users/myles/R/Deconv/brain_atlas.rdata",
#      list = c("mk", "mkNN", "mkm"))

# reload
load("/Users/myles/R/Deconv/brain_atlas.rdata")

# all except
gp <- colnames(mk$groupmeans)
except <- setdiff  # use with pipe

mks <- updateMarkers(mk, remove_group = gp |> except("Splatter"))

# simulate bulk
mk <- updateMarkers(mk, expfilter = 0.2, nsubclass = 20)

set.seed(3)
sim_counts <- generate_samples(mk, 30)
sim_percent <- sim_counts / rowSums(sim_counts) * 100

# tougher
sim_sampled <- simulate_bulk(mat, sim_counts, meta$roi,
                             times = 1, method = "dirichlet")
# 35 mins (Intel)
rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)

#############
# set up bulk simulation for non-neuronal
set.seed(3)
sim_countsNN <- generate_samples(mkNN, 30)
sim_percentNN <- sim_countsNN / rowSums(sim_countsNN) * 100

# tougher
sim_sampledNN <- simulate_bulk(mat2, sim_countsNN, meta2$cell_type,
                             times = 1, method = "dirichlet")
# 6.32 mins (Intel)
rownames(sim_sampledNN) <- gene2symbol(rownames(sim_sampledNN), ensDb_v110)

# merge
identical(rownames(sim_sampled), rownames(sim_sampledNN))

sim_sampled_merge <- sim_sampled + sim_sampledNN
sim_percent_merge <- sim_counts_merge / rowSums(sim_counts_merge) * 100
sim_counts_merge <- cbind(sim_counts, sim_countsNN)

mk$cell_table
mkNN$cell_table
mkm$cell_table

mkm <- updateMarkers(mkm, expfilter = 0.2, nsubclass = 20)

fitm <- deconvolute(mkm, sim_sampled_merge,
                    weight_method = "equal",
                    arith_mean = T,
                    use_filter = F, cores = 8)

mset <- metric_set(sim_percent_merge, fitm$subclass$percent)
summary(mset)
summary(fitm$subclass$comp_amount)

pdf("/users/myles/dropbox/R scripts/Deconv/pdf/brain_sim_merge_dirich_eq.pdf",
    width = 12, height = 12.5)
plot_set(sim_counts_merge, fitm$subclass$output, show_zero = T,
         mfrow = c(8, 8))
plot_set(sim_percent_merge, fitm$subclass$percent, show_zero = T,
         mfrow = c(8, 8))
dev.off()

save(file = "/Users/myles/R/Deconv/brain_atlas2.rdata",
     list = c("mk", "mkNN", "mkm",
              "sim_counts", "sim_percent", "sim_sampled",
              "sim_countsNN", "sim_percentNN", "sim_sampledNN",
              "sim_counts_merge", "sim_percent_merge", "sim_sampled_merge"))

# reload
load("/Users/myles/R/Deconv/brain_atlas2.rdata")

qqm <- quantile_map(mk, mkNN)
pdf("/Users/myles/dropbox/R scripts/Deconv/pdf/brain_qq.pdf",
    width = 5, height = 5)
plot(qqm, tcl = -0.3, mgp = c(1.6,0.5,0), cex = 0.7,
     xlab = "Neurons", ylab = "Non-neuronal cells", bty = "l")
dev.off()

## start here
noise_check <- function(noiseFUN, noise_set, reps = 3) {
  res <- lapply(noise_set, function(noi) {
    out <- vapply(1:reps, function(j) {
      sim_noise <- do.call(noiseFUN, list(counts = sim_sampled_merge, sd = noi))
      fit <- deconvolute(mkm, sim_noise,
                         arith_mean = T, use_filter = F, comp_amount = 1, cores = 16)
      err <- fit$subclass$output - sim_counts_merge
      rmse <- sqrt(colMeans(err^2))
      vapply(1:6, function(i) {
        se <- fit$subclass[[paste0("se", i)]]
        mse <- if (is.matrix(se)) sqrt(colMeans(se^2)) else se
        Rsq(rmse, mse)
      }, numeric(1))
    }, numeric(6))
    matrix(c(rowMeans(out), matrixStats::rowSds(out) / sqrt(reps)), ncol = 2,
           dimnames = list(NULL, c("mean", "SEM")))
  })
  
  res <- do.call(rbind, res)
  res2 <- data.frame(res)
  res2$method <- factor(c("OLS", "OLS.v2", "HC0", "HC2", "HC3", "var.e"),
                        levels = c("OLS", "OLS.v2", "HC0", "HC2", "HC3", "var.e"))
  res2$noise <- rep(noise_set, each = 6)
  res2
}

plot_noise <- function(res2, main = "") {
  p <- ggplot(res2, aes(x = noise, y = mean, colour = method)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM),
                  width = diff(range(res2$noise)) *0.02) +
    ggtitle(main) +
    ylab(expression(R^2)) +
    ylim(NA, 1) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          plot.title = element_text(size = 10))
  p
}


r1 <- noise_check("add_noise", noise_set = 0:5 * 1000)
r2 <- noise_check("log_noise", noise_set = 0:5 / 20)
r3 <- noise_check("sqrt_noise", noise_set = 0:5 * 1000)
r4 <- noise_check("shift_noise", noise_set = 0:5 / 10)

p1 <- plot_noise(r1, main = "Gaussian noise")
p2 <- plot_noise(r2, main = "log noise")
p3 <- plot_noise(r3, main = "sqrt noise")
p4 <- plot_noise(r4, main = "shift noise")


library(ggpubr)
pdf("/users/myles/dropbox/R scripts/deconv/pdf/brain_compare_SE_noise.pdf",
    width = 8.5)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE,
          legend = "right", align = "hv")
dev.off()

pdf("/users/myles/dropbox/R scripts/deconv/pdf/brain_compare_SE_noise horiz.pdf",
    width = 12.5, height = 3)
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = TRUE,
          legend = "right", align = "hv")
dev.off()

pdf("/users/myles/dropbox/R scripts/deconv/pdf/brain_compare_SE_noise horiz 0_1 only.pdf",
    width = 12.5, height = 3)
ggarrange(p1 +ylim(0,1), p2 +ylim(0,1), p3 +ylim(-6,1), p4 +ylim(0,1),
          ncol = 4, nrow = 1, common.legend = TRUE,
          legend = "right", align = "hv")
dev.off()

err <- fitm$subclass$output - sim_counts_merge
rmse <- sqrt(colMeans(err^2))
mse <- colMeans(fitm$subclass$se4)
mse <- fitm$subclass$se6

plot(rmse, mse, las = 1, bty = "l", tcl = -0.3, mgp = c(2, 0.5, 0))
abline(0, 1, col = "red")
abline(0, 1/1.96, col = "purple")
Rsq(rmse, mse)
sort(mse)
fitm$subclass$resvar

pdf("/users/myles/dropbox/R scripts/deconv/pdf/brain_SE_v_RMSE.pdf",
    width = 6.5, height = 6)
op <- par(mar = c(4, 4, 2, 1))
plot(rmse, mse, las = 1, bty = "l", tcl = -0.3, mgp = c(2.5, 0.5, 0),
     xlab = "Actual error (RMSE)", ylab = "Deconvolution std error")
abline(0, 1, col = "red")
par(op)
dev.off()

save(file = "/Users/myles/R/Deconv/brain_SE_noise.rdata",
     list = c("mk", "mkNN", "mkm",
              "sim_counts", "sim_percent", "sim_sampled",
              "sim_countsNN", "sim_percentNN", "sim_sampledNN",
              "sim_counts_merge", "sim_percent_merge", "sim_sampled_merge",
              "r1", "r2", "r3", "r4",
              "p1", "p2", "p3", "p4"))
