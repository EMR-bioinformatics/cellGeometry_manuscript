# Tabula sapiens full dataset

# https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5

# full dataset
# 45 GB
# 1,136,218 cells

library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)

# Create SingleCellExperiment object that points to on-disk h5ad file

library(devtools)
setwd("/users/myles/documents/github/cellGeometry")
load_all()

library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]

# try full dataset
tfull <- readH5AD("/Users/myles/R/Deconv/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad",
                  use_hdf5 = TRUE, reader = "R")

mat <- tfull@assays@data$X
rownames(mat) <- rownames(tfull)  # need to add rownames (genes)
meta <- tfull@colData@listData

mf <- cellMarkers(mat, subclass = meta$cell_type, 
                  cellgroup = meta$broad_cell_class,
                  dual_mean = T, cores = 7)

# using DelayedArray::rowMeans
# 6.63 + 6.31 mins (intel)
# 4.64 + 3.8 mins (intel)
# 1.8 + 2.14 mins (M3)

mf <- gene2symbol(mf, ensDb_v110)

saveRDS(mf, "/Users/myles/R/Deconv/tabula_full_mk.rds")
mk <- readRDS("/Users/myles/R/Deconv/tabula_full_mk.rds")

set.seed(3)
sim_counts <- generate_samples(mk, 50)
sim_percent <- sim_counts / rowSums(sim_counts) * 100

# sampled dirichlet
# 8.6 mins on M3
set.seed(99)
sim_sampled <- simulate_bulk(mat, sim_counts, meta$cell_type, times = 1,
                             method = "dirichlet")
rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)


fit2 <- deconvolute(mk, sim_sampled,
                    arith_mean = T, use_filter = F, cores = 12)

plot_set(sim_counts, fit2$subclass$output, show_zero = T, show_identity = T)
mset <- metric_set(sim_percent, fit2$subclass$percent)
summary(mset)

err <- fit2$subclass$output - sim_counts
rmse <- sqrt(colMeans(err^2))
mse <- colMeans(fit2$subclass$se4)
mse <- fit2$subclass$se

plot(rmse, mse, las = 1, bty = "l", tcl = -0.3, mgp = c(2, 0.5, 0))
abline(0, 1, col = "red")
abline(0, 1/1.96, col = "purple", lty = 2)

dat <- data.frame(rmse, mse)
library(easylabel)
easylabel(dat, "rmse", "mse", zeroline = F,
          colScheme = "blue", alpha = 0.7,
          panel.last = {abline(0, 1, col = "red")})

plot_se <- function(fit, sim_counts) {
  err <- fit$subclass$output - sim_counts
  rmse <- sqrt(colMeans(err^2))
  mse <- fit2$subclass$se6
  mse <- sqrt(colMeans(fit2$subclass$se5^2))
  plot(rmse, mse, las = 1, bty = "l", tcl = -0.3, mgp = c(2, 0.5, 0))
  abline(0, 1, col = "red")
  # abline(0, 0.5, col = "purple")
  # lmfit <- lm(mse ~ rmse)
  # print(lmfit)
}
plot_se(fit2, sim_counts)

noise_check <- function(noiseFUN, noise_set, reps = 3) {
  res <- lapply(noise_set, function(noi) {
    out <- vapply(1:reps, function(j) {
      sim_noise <- do.call(noiseFUN, list(counts = sim_sampled, sd = noi))
      fit2 <- deconvolute(mk, sim_noise,
                          arith_mean = T, use_filter = F, comp_amount = 1, cores = 16)
      err <- fit2$subclass$output - sim_counts
      rmse <- sqrt(colMeans(err^2))
      vapply(1:6, function(i) {
        se <- fit2$subclass[[paste0("se", i)]]
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
                        c("OLS", "OLS.v2", "HC0", "HC2", "HC3", "var.e"))
  res2$noise <- rep(noise_set, each = 6)
  list(res = res2, noise_set = noise_set)
}

noise_plot <- function(obj, main = "") {
  res <- obj$res
  noise_set <- obj$noise_set
  p <- ggplot(res, aes(x = noise, y = mean, colour = method)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM),
                  width = diff(range(noise_set)) *0.02) +
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

p1 <- noise_plot(r1, main = "Gaussian noise")
p2 <- noise_plot(r2, main = "log noise")
p3 <- noise_plot(r3, main = "sqrt noise")
p4 <- noise_plot(r4, main = "shift noise")

# show var.e only
p1 + ylim(-0.5, 1)

library(ggpubr)
pdf("/users/myles/dropbox/R scripts/deconv/pdf/tabula_compare_SE_noise.pdf",
    width = 8.5)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE,
          legend = "right", align = "hv")
dev.off()

pdf("/users/myles/dropbox/R scripts/deconv/pdf/tabula_compare_SE_noise horiz.pdf",
    width = 12.5, height = 3)
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = TRUE,
          legend = "right", align = "hv")
dev.off()

pdf("/users/myles/dropbox/R scripts/deconv/pdf/tabula_compare_SE_noise horiz 0_1 only.pdf",
    width = 12.5, height = 3)
ggarrange(p1 + ylim(-1.2, 1), p2 + ylim(-1.2, 1), p3 + ylim(-1.2, 1), p4 + ylim(-1.2, 1), ncol = 4, nrow = 1, common.legend = TRUE,
          legend = "right", align = "hv")
dev.off()


# plot heteroscedasticity
gm <- sim_sampled[rownames(fit2$subclass$residuals), ]
plot(gm, fit2$subclass$residuals, log = "x",
     cex = 0.7, col = adjustcolor("black", 0.2), pch = 16, bty = "l")
abline(1,0, col = "red")

plot(rowMeans(gm), matrixStats::rowVars(fit2$subclass$residuals),
     log = "x", cex = 0.7, col = adjustcolor("black", 0.2), pch = 16, bty = "l")

plotly::plot_ly(x = rowMeans(gm), y = matrixStats::rowVars(fit2$subclass$residuals),
                type = "scattergl", text = rownames(gm))

hist(scale(log2(fit2$subclass$var.e))[,1])
sort(scale(log2(fit2$subclass$var.e))[,1], decreasing = T)[1:10]
hist(rstudent(fit2))


save(file = "/Users/myles/R/Deconv/tabula_sim.rdata",
     list = c("mk", "sim_counts", "sim_percent", "sim_sampled"))

save(file = "/Users/myles/R/Deconv/tabula_SE_noise.rdata",
     list = c("mk", "sim_counts", "sim_percent", "sim_sampled",
              "r1", "r2", "r3", "r4",
              "p1", "p2", "p3", "p4"))

