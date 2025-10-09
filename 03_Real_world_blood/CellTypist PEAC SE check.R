# Cell Typist

# https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3

library(devtools)
library(AnnotationHub)

setwd("/users/myles/documents/github/cellGeometry")
load_all()

ah <- AnnotationHub()
# ah <- AnnotationHub(localHub = T)
ensDb_v110 <- ah[["AH113665"]]

# h5ad
library(zellkonverter)
library(SingleCellExperiment)

typist_h5 <- readH5AD("/Users/myles/R/Deconv/2ac906a5-9725-4258-8e36-21a9f6c0302a.h5ad",
                    use_hdf5 = TRUE, reader = "R")

mat <- typist_h5@assays@data$X
rownames(mat) <- rownames(typist_h5)
meta <- typist_h5@colData@listData
# table(meta$Majority_voting_CellTypist)

subcl <- meta$Majority_voting_CellTypist
cellgrp <- meta$Majority_voting_CellTypist_high

subcl[meta$tissue != "blood"] <- NA
cellgrp[meta$tissue != "blood"] <- NA

mk <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp,
                  dual_mean = T, cores = 8)
mk <- gene2symbol(mk, ensDb_v110)

mk <- updateMarkers(mk,
                    remove_subclass = c("Helper T cells", "Cytotoxic T cells"))

# full 43 subclasses
subcl <- meta$Majority_voting_CellTypist
cellgrp <- meta$Majority_voting_CellTypist_high
mk0 <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp,
                   dual_mean = T, cores = 8)
mk0 <- gene2symbol(mk0, ensDb_v110)
mk0 <- updateMarkers(mk0,
                     remove_subclass = c("Helper T cells", "Cytotoxic T cells"))

# simulation
set.seed(3)
sim_counts <- generate_samples(mk, 50)
sim_percent <- sim_counts / rowSums(sim_counts) * 100

# sampled dirichlet
set.seed(99)
sim_sampled <- simulate_bulk(mat, sim_counts, subcl, times = 1,
                             method = "dirichlet")
rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)

fit0 <- deconvolute(mk, sim_sampled, arith_mean = T, use_filter = F, cores = 8, npass = 2)
hist(scale(log2(fit0$subclass$var.e +1)))

barplot(sort(fit0$subclass$se), horiz = T, las = 1, cex.names = 0.5)
head(sort(fit0$subclass$hat, decreasing = T), 10)

# real bulk blood data
PEAC.blood <- readRDS("/users/myles/R/deconv/peac_bld_counts.rds")
mk <- updateMarkers(mk, bulkdata = PEAC.blood)

mk <- updateMarkers(mk, bulkdata = PEAC.blood, nsubclass = 30)
mk <- updateMarkers(mk, bulkdata = PEAC.blood, nsubclass = 25)
mk <- updateMarkers(mk, bulkdata = PEAC.blood, nsubclass = 5000,
                    expfilter = 0.1)
mk <- updateMarkers(mk, bulkdata = PEAC.blood, nsubclass = 500,
                    expfilter = 0.2)

# from nsubclass 29 to 30
barplot(sort(mk$genemeans["HBB",]), horiz = T, las = 1)

fit <- deconvolute(mk, PEAC.blood, cores = 8, npass = 1)
fit <- deconvolute(mk, PEAC.blood, cores = 8, npass = 2)
fit <- deconvolute(mk, PEAC.blood, cores = 8, npass = 5)

# test multipass options
fit <- deconvolute(mk, PEAC.blood, cores = 8, npass = 5, outlier_cutoff = 3)

fit <- deconvolute(mk, PEAC.blood, cores = 8, npass = 3,
                   outlier_method = "cooks", outlier_cutoff = 0.5,
                   outlier_quantile = 0.9)
fit <- deconvolute(mk, PEAC.blood, cores = 8, npass = 5,
                   outlier_method = "rstud", outlier_cutoff = 4,
                   outlier_quantile = 0.9)

barplot(sqrt(sort(fit$subclass$var.e, decreasing = T)[1:10]))

fit0 <- deconvolute(mk, sim_sampled, arith_mean = T, use_filter = F, cores = 8,
                    npass = 1)
fit0 <- deconvolute(mk, sim_sampled, arith_mean = T, use_filter = F, cores = 8,
                    npass = 5, outlier_method = "rstud")

hist(fit$subclass$se)
barplot(fit$subclass$comp_amount, las = 3, cex.names = 0.7)

hist(fit$subclass$var.e, 100)
hist(log2(fit$subclass$var.e +1))
hist(sqrt(fit$subclass$var.e), 100)
hist(scale(log2(fit$subclass$var.e +1)))
range(scale(log2(fit$subclass$var.e +1)))
sort(scale(log2(fit$subclass$var.e +1))[,1], decreasing = TRUE)[1:10]

# HBG2    ACSL1   NBEAL2    LRRK2     FPR1    MEGF9    ANPEP   IFI44L     GBP5   CSF2RB 
# 3.296938 3.223698 3.199919 3.090549 2.936006 2.925309 2.845193 2.831031 2.814257 2.809591 
# ACSL1, LRRK2, FPR1, MGEF9 highly expressed in neutrophils

# plot actual vs predicted gene expression
cellmat <- 2^fit$mk$genemeans -1
output <- fit$subclass$output
pred_exp <- tcrossprod(cellmat, output)
overlap_genes <- intersect(rownames(pred_exp), rownames(PEAC.blood))

plot(log2(PEAC.blood[overlap_genes, ] +1), log2(pred_exp[overlap_genes, ] +1),
     cex = 0.3, col = adjustcolor("black", 0.04), pch = 16,
     las = 1, xlab = "real bulk", ylab = "inferred", bty = "l",
     tcl = -0.4, mgp = c(2, 0.6, 0))
abline(0,1, col = "red")

smoothScatter(log2(PEAC.blood[overlap_genes, ] +1), log2(pred_exp[overlap_genes, ] +1),
     las = 1, xlab = "real bulk", ylab = "inferred", bty = "l",
     tcl = -0.4, mgp = c(2, 0.6, 0))

rstud <- rstandard(fit0)
rstud <- rstudent(fit0)
sort(rowMaxs(abs(rstud)), decreasing = T)[1:10]
sort(rowMeans(abs(rstud)), decreasing = T)[1:10]

plot(rowMeans(abs(rstud)), fit0$subclass$var.e)
plot(scale(log2(fit0$subclass$var.e)), rowMeans(abs(rstud)))
# simulation

mk <- updateMarkers(mk, nsubclass = 25)
fit0 <- deconvolute(mk, sim_sampled, arith_mean = T, use_filter = F, cores = 8, npass = 2)

pdf("/users/myles/dropbox/R scripts/deconv/pdf/typist sim var_e vs rstud or cooks.pdf",
    width = 9.5, height = 4.5)
op <- par(mar = c(4, 4, 2, 1), mfrow= c(1, 2))
plot(sqrt(fit0$subclass$var.e), rowMeans(abs(rstudent(fit0))),
     cex = 0.8, col = adjustcolor("black", 0.5), pch = 16,
     las = 1, xlab = "sd[e]", ylab = "mean abs Studentized residuals", bty = "l",
     tcl = -0.4, mgp = c(2, 0.6, 0))
abline(lm(rowMeans(abs(rstud)) ~ sqrt(fit0$subclass$var.e)),
       col = "red", lwd = 1.5)
cor1 <- cor.test(sqrt(fit0$subclass$var.e), rowMeans(abs(rstudent(fit0))), method = "spear", exact = F)
mtext(bquote(rho ~"=" ~ .(format(cor1$estimate, digits = 3)) ~ ", P =" ~ .(format(cor1$p.value, digits = 3))), adj = 0)

plot(sqrt(fit0$subclass$var.e), rowMeans(cooks.distance(fit0)),
     cex = 0.8, col = adjustcolor("black", 0.5), pch = 16,
     xlab = "", ylab = "mean Cook's distance", bty = "l",
     tcl = -0.3, mgp = c(2.9, 0.4, 0),
     las = 1)
mtext("sd[e]", side = 1, line = 2)
lmf <- lm(rowMeans(cooks.distance(fit0)) ~ sqrt(fit0$subclass$var.e))
abline(lmf, col = "red", lwd = 1.5)
cor2 <- cor.test(sqrt(fit0$subclass$var.e), rowMeans(cooks.distance(fit0)), method = "spearman", exact = F)
mtext(bquote(rho ~"=" ~ .(format(cor2$estimate, digits = 3)) ~ ", P =" ~ .(format(cor2$p.value, digits = 3))), adj = 0)
par(op)
dev.off()


# try same with PEAC blood
# var.e vs rstud or cooks
pdf("/users/myles/dropbox/R scripts/deconv/pdf/typist PEAC_blood fit5 var_e vs rstud or cooks.pdf",
    width = 9.5, height = 4.5)
op <- par(mar = c(4, 4, 2, 1), mfrow= c(1, 2))
plot(sqrt(fit5$subclass$var.e), rowMeans(abs(rstudent(fit5))),
     cex = 0.8, col = adjustcolor("black", 0.5), pch = 16,
     las = 1, xlab = "sd[e]", ylab = "mean abs Studentized residuals", bty = "l",
     tcl = -0.4, mgp = c(2, 0.6, 0))
abline(lm(rowMeans(abs(rstudent(fit5))) ~ sqrt(fit5$subclass$var.e)),
       col = "red", lwd = 1.5)
cor1 <- cor.test(sqrt(fit5$subclass$var.e), rowMeans(abs(rstudent(fit5))), method = "spear", exact = F)
mtext(bquote(rho ~"=" ~ .(format(cor1$estimate, digits = 3)) ~ ", P =" ~ .(format(cor1$p.value, digits = 3))), adj = 0)

plot(sqrt(fit5$subclass$var.e), rowMeans(cooks.distance(fit5)),
     cex = 0.8, col = adjustcolor("black", 0.5), pch = 16,
     xlab = "", ylab = "mean Cook's distance", bty = "l",
     tcl = -0.3, mgp = c(2.9, 0.4, 0),
     las = 1)
mtext("sd[e]", side = 1, line = 2)
lmf <- lm(rowMeans(cooks.distance(fit5)) ~ sqrt(fit5$subclass$var.e))
abline(lmf, col = "red", lwd = 1.5)
cor2 <- cor.test(sqrt(fit5$subclass$var.e), rowMeans(cooks.distance(fit5)), method = "spearman", exact = F)
mtext(bquote(rho ~"=" ~ .(format(cor2$estimate, digits = 3)) ~ ", P =" ~ .(format(cor2$p.value, digits = 3))), adj = 0)
par(op)
dev.off()


# note
# sd[e] is sd(weighted residuals) per gene


# actual blood
plot(sqrt(fit$subclass$var.e), rowMeans(abs(rstudent(fit))))
sort(fit$subclass$var.e, decreasing = T)[1:10]
sort(scale(log2(fit$subclass$var.e))[,1], decreasing = T)[1:20]
# rstudent cutoff is 10
# var.e cutoff is 4 (on Z score of log2 of var.e)

plotly::plot_ly(x = sqrt(fit$subclass$var.e), y = rowMeans(abs(rstudent(fit))),
                text = names(fit$subclass$var.e),
                type = "scattergl")

# cooks
plot(sqrt(fit$subclass$var.e), rowMeans(cooks.distance(fit)))

plotly::plot_ly(x = sqrt(fit$subclass$var.e), y = rowMeans(cooks.distance(fit)),
                text = names(fit$subclass$var.e),
                type = "scattergl")

plotly::plot_ly(x = rowMeans(abs(rstudent(fit))), y = rowMeans(cooks.distance(fit)),
                text = names(fit$subclass$var.e),
                type = "scattergl")

ckd <- cooks.distance(fit0)
range(ckd)
sort(rowMaxs(ckd), decreasing = T)[1:10]
sort(rowMeans(ckd), decreasing = T)[1:10]

# OLS comparison
fit_ols <- ols(mk, PEAC.blood)
fit_ols <- ols(mk, sim_sampled)
fit_ols <- ols(mk, sim_sampled, arith_mean = T, use_filter = F)
fit_ols$fit

hist(residuals(fit_ols$fit))
hist(residuals(fit_ols$fitw))

sort(matrixStats::rowQuantiles(rstud, probs = 0.9), decreasing = T)[1:10]
sort(matrixStats::rowQuantiles(ckd, probs = 0.9), decreasing = T)[1:10]

ols_rstud <- rstudent(fit_ols$fit)
ols_rstud <- rstudent(fit_ols$fitw)
range(ols_rstud)
hist(ols_rstud, 100)
sort(rowMaxs(abs(ols_rstud)), decreasing = T)[1:10]
sort(rowMeans(abs(ols_rstud)), decreasing = T)[1:10]
sort(matrixStats::rowQuantiles(ols_rstud, probs = 0.9), decreasing = T)[1:10]

ols_cooks <- cooks.distance(fit_ols$fit)
ols_cooks <- cooks.distance(fit_ols$fitw)
range(ols_cooks)
sort(rowMaxs(ols_cooks), decreasing = T)[1:10]
sort(rowMeans(ols_cooks), decreasing = T)[1:10]

ols_dffits <- dffits(fit_ols$fit)
range(ols_dffits)
sort(rowMaxs(ols_dffits), decreasing = T)[1:10]
sort(rowMeans(ols_dffits), decreasing = T)[1:10]

# noise test, full dataset
set.seed(3)
sim_counts <- generate_samples(mk0, 50)
sim_percent <- sim_counts / rowSums(sim_counts) * 100

# sampled dirichlet
set.seed(99)
sim_sampled <- simulate_bulk(mat, sim_counts, subcl, times = 1,
                             method = "dirichlet")
rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)

sim_noise <- sim_sampled
sim_noise <- add_noise(sim_sampled, sd = 1000)
sim_noise <- sqrt_noise(sim_sampled, sd = 1000)
sim_noise <- log_noise(sim_sampled, sd = 0.5)
sim_noise <- shift_noise(sim_sampled, sd = 0.1)
sim_noise <- shift_noise(sim_sampled, sd = 2, p = 0.3)

# plot noise
plot(sim_sampled, sim_noise, log = "xy", las = 1,
     cex = 0.6, col = adjustcolor("black", 0.12), pch = 16, bty = "l")

mk0 <- updateMarkers(mk0, nsubclass = 200)

fit2 <- deconvolute(mk0, sim_noise,
                    arith_mean = T, use_filter = F, cores = 8,
                    npass = 1, outlier_method = "rstud",
                    outlier_cutoff = 2.5)

plot_set(sim_counts, fit2$subclass$output)
summary(metric_set(sim_counts, fit2$subclass$output))
hist(fit2$subclass$se)

shift_set <- attr(sim_noise, "shift")[rownames(sim_noise) %in% mk0$geneset]
names(shift_set) <- mk0$geneset
nonzero_shift <- names(shift_set)[shift_set != 0]
intersect(nonzero_shift, names(fit2$subclass$removed))

res <- tune_deconv(mk0, sim_noise, sim_counts,
                   grid = list(nsubclass = c(25, 50, 100, 200),
                               npass = 1:5,
                               outlier_method = c("var.e", "cook", "rstud")),
                   arith_mean = T, use_filter = F,
                   cores = 16)

plot_tune(res, xvar = "nsubclass", group = "npass")
plot_tune(res, xvar = "nsubclass", group = "outlier_method")
plot_tune(res, xvar = "nsubclass", group = "outlier_method", fix = list(npass = 2))
plot_tune(res, xvar = "npass", group = "outlier_method")
plot_tune(res, xvar = "npass", group = "outlier_method", fix = list(nsubclass = 200))

res <- tune_deconv(mk0, sim_noise, sim_counts,
                   grid = list(# nsubclass = c(25, 50, 100, 200),
                               npass = c(1:5),
                               # outlier_cutoff = c(2, 4, 6, 8, 10)
                               outlier_cutoff = c(2, 2.5, 3, 3.5, 4)
                               # outlier_quantile = c(0.5, 0.9, 1)
                               ),
                   arith_mean = T, use_filter = F, outlier_method = "rstud",
                   cores = 16)

plot_tune(res, xvar = "npass", group = NULL)
plot_tune(res, xvar = "npass", group = "outlier_cutoff")
plot_tune(res, xvar = "outlier_cutoff", group = "outlier_quantile")
plot_tune(res, xvar = "outlier_cutoff", group = "npass")
plot_tune(res, xvar = "ngene", group = NULL)

unique(res$ngene)

# plot SE vs RMSE
err <- fit0$subclass$output - sim_counts
rmse <- sqrt(colMeans(err^2))
# mse <- colMeans(fit2$subclass$se4)
mse <- fit0$subclass$se

str(fit2$subclass)

plot(rmse, mse, las = 1, bty = "l", tcl = -0.3, mgp = c(2, 0.5, 0))
abline(0, 1, col = "red")

mk <- updateMarkers(mk, bulkdata = PEAC.blood, nsubclass = 500,
                    expfilter = 0.2)

# plot residuals & heteroscedasticity
fit <- deconvolute(mk, PEAC.blood, cores = 8, npass = 1)
fit2 <- deconvolute(mk, PEAC.blood, cores = 8, npass = 2)
fit3 <- deconvolute(mk, PEAC.blood, cores = 8, npass = 3)
fit5 <- deconvolute(mk, PEAC.blood, cores = 8, npass = 5)

pdf("/users/myles/dropbox/R scripts/deconv/pdf/typist sim resid.pdf",
    width = 10, height = 3)
op <- par(mfrow = c(1, 3), mgp = c(2.5, 0.4, 0), tcl = -0.3, las = 1)
plot_residuals(fit0, sim_sampled)
plot_residuals(fit0, sim_sampled, type = "weight")
plot_residuals(fit0, sim_sampled, type = "student")
par(op)
dev.off()

# full plot for supplement
pdf("/users/myles/dropbox/R scripts/deconv/pdf/typist PEAC blood resid vs npass.pdf",
    width = 10, height = 12)
op <- par(mar = c(4, 5, 1.5, 1.1), mfrow = c(4, 3), mgp = c(3.5, 0.4, 0),
          tcl = -0.3, las = 1)
plot_residuals(fit, PEAC.blood)
plot_residuals(fit, PEAC.blood, type = "weight")
plot_residuals(fit, PEAC.blood, type = "student")
plot_residuals(fit2, PEAC.blood)
plot_residuals(fit2, PEAC.blood, type = "weight")
plot_residuals(fit2, PEAC.blood, type = "student")
plot_residuals(fit3, PEAC.blood)
plot_residuals(fit3, PEAC.blood, type = "weight")
plot_residuals(fit3, PEAC.blood, type = "student")
plot_residuals(fit5, PEAC.blood)
plot_residuals(fit5, PEAC.blood, type = "weight")
plot_residuals(fit5, PEAC.blood, type = "student")
par(op)
dev.off()

# studentized only
pdf("/users/myles/dropbox/R scripts/deconv/pdf/typist PEAC blood rstud vs npass.pdf",
    width = 11, height = 2.6)
op <- par(mar = c(3.5, 3.5, 1.5, 1), mfrow = c(1, 4), mgp = c(2, 0.4, 0),
          tcl = -0.3, las = 1)
plot_residuals(fit, PEAC.blood, type = "student")
mtext("Pass 1", adj = 0.02, cex = par("cex"))
plot_residuals(fit2, PEAC.blood, type = "student")
mtext("Pass 2", adj = 0.02, cex = par("cex"))
plot_residuals(fit3, PEAC.blood, type = "student")
mtext("Pass 3", adj = 0.02, cex = par("cex"))
plot_residuals(fit5, PEAC.blood, type = "student", xlab = "a")
mtext("Pass 4", adj = 0.02, cex = par("cex"))
par(op)
dev.off()

ggplot_residuals(fit, PEAC.blood, type = "student")

save.image("/users/myles/R/deconv/typist_PEAC_SE_check.rdata")
