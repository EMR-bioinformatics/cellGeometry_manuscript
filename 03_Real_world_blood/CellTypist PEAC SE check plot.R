####Background####

# From Myles' script 'CellTypist PEAC SE check.R'
# To obtain the 'typist sim resid' plot and the 'typist PEAC blood rstud vs npass' plot

####simulation####

library(ggplot2)
library(ggrastr)

ggplot_residualsv2 <- function(fit, test, type = c("reg", "student", "weight"),
                               show_outliers = TRUE) {
  nm <- fit$call$test
  .call <- match.call()
  if (nm != .call$test) {
    message("`", .call$test, "` does not match the bulk dataset `", nm,
            "` used in `", .call$fit, "`")
  }
  dat <- plot_residuals(fit, test, type, show_outliers, show_plot = FALSE)
  type <- match.arg(type)
  dat$expr[dat$expr < 1] <- NA
  ylab <- switch(type, student = "Studentized residuals",
                 weight = "Weighted residuals",
                 "Raw residuals")
  ggplot(dat, aes(x = .data$expr, y = .data$res)) +
    rasterize(geom_point(aes(colour = .data$outlier), na.rm = TRUE),
              #dpi = 1000,
              dpi = 400,
              scale = 0.6) +
    geom_hline(yintercept = 0, color = "blue") +
    scale_x_log10(guide = "axis_logticks") +
    scale_colour_manual(values = c(adjustcolor("black", 0.2),
                                   adjustcolor("red", 0.7))) +
    xlab("Bulk gene expression") + ylab(ylab) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.title = element_text(size = 11))
}

library(ggpubr)

#load with Myles' data

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/Myles_script/typist_PEAC_SE_check.rdata")

pdf("typist sim resid_Myles.pdf",
    width = 10, height = 3)
ggarrange(
  ggplot_residualsv2(fit0, sim_sampled),
  ggplot_residualsv2(fit0, sim_sampled, type = "weight"), 
  ggplot_residualsv2(fit0, sim_sampled, type = "student"),
  ncol = 3)
dev.off()

pdf("typist sim var_e vs rstud or cooks Myles.pdf",
    width = 7.6, height = 3.6)
op <- par(mar = c(4, 4, 2, 1), mfrow= c(1, 2))
plot(sqrt(fit0$subclass$var.e), rowMeans(abs(rstudent(fit0))),
     #cex = 0.8,
     cex = 0.64, col = adjustcolor("black", 0.5), pch = 16,
     las = 1, xlab = "sd[e]", ylab = "mean abs Studentized residuals", bty = "l",
     tcl = -0.4, mgp = c(2, 0.6, 0))
abline(lm(rowMeans(abs(rstud)) ~ sqrt(fit0$subclass$var.e)),
       col = "red", lwd = 1.5)
cor1 <- cor.test(sqrt(fit0$subclass$var.e), rowMeans(abs(rstudent(fit0))), method = "spear", exact = F)
mtext(bquote(rho ~"=" ~ .(format(cor1$estimate, digits = 3)) ~ ", P =" ~ .(format(cor1$p.value, digits = 3))), adj = 0)

plot(sqrt(fit0$subclass$var.e), rowMeans(cooks.distance(fit0)),
     #cex = 0.8,
     cex = 0.64, col = adjustcolor("black", 0.5), pch = 16,
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

####real world PEAC blood####

library(cellGeometry)

PEAC.blood <- readRDS("peac_bld_counts.rds")

mk_PEAC <- updateMarkers(mk, bulkdata = PEAC.blood, nsubclass = 500,
                         expfilter = 0.2)

fit <- deconvolute(mk_PEAC, PEAC.blood, cores = 8, npass = 1)
fit2 <- deconvolute(mk_PEAC, PEAC.blood, cores = 8, npass = 2)
fit3 <- deconvolute(mk_PEAC, PEAC.blood, cores = 8, npass = 3)
fit5 <- deconvolute(mk_PEAC, PEAC.blood, cores = 8, npass = 5)

#stopped at npass = 4
#HBA1, HBB 
#IFIT2, IFIT3, OAS3, VNN2, AQP9, LTF, FFAR2, DYSF, PTGS2, FPR2 
#DEFA3 

# save(mk_PEAC, fit, fit2, fit3, fit5,
#      file = "typist_PEAC_bld.rdata")

pdf("typist PEAC blood rstud vs npass.pdf",
    width = 11, height = 2.6)
ggarrange(
ggplot_residualsv2(fit, PEAC.blood, type = "student") + ggtitle("Pass = 1"),
ggplot_residualsv2(fit2, PEAC.blood, type = "student") + ggtitle("Pass = 2"),
ggplot_residualsv2(fit3, PEAC.blood, type = "student") + ggtitle("Pass = 3"),
ggplot_residualsv2(fit5, PEAC.blood, type = "student") + ggtitle("Pass = 4"),
ncol = 4)
dev.off()

#for supplement
pdf("typist PEAC blood resid vs npass supp v2.pdf",
    width = 10, height = 12)
ggarrange(
  ggplot_residualsv2(fit, PEAC.blood),
  ggplot_residualsv2(fit, PEAC.blood, type = "weight"),
  ggplot_residualsv2(fit, PEAC.blood, type = "student"),
  ggplot_residualsv2(fit2, PEAC.blood),
  ggplot_residualsv2(fit2, PEAC.blood, type = "weight"),
  ggplot_residualsv2(fit2, PEAC.blood, type = "student"),
  ggplot_residualsv2(fit3, PEAC.blood),
  ggplot_residualsv2(fit3, PEAC.blood, type = "weight"),
  ggplot_residualsv2(fit3, PEAC.blood, type = "student"),
  ggplot_residualsv2(fit5, PEAC.blood),
  ggplot_residualsv2(fit5, PEAC.blood, type = "weight"),
  ggplot_residualsv2(fit5, PEAC.blood, type = "student"),
  ncol = 3, nrow = 4)
dev.off()


####SE with npass####

all(names(fit5$subclass$se) == names(fit2$subclass$se))
#TRUE

SE <- data.frame("SE" = c(fit$subclass$se,
                          fit2$subclass$se,
                          fit3$subclass$se,
                          fit5$subclass$se),
                 "subclass" = names(fit$subclass$se),
                 "npass" = rep(c(1, 2, 3, 4),
                               each = length(fit$subclass$se)))

col_scheme <- c("Memory B" = "skyblue1",
                "Naive B" = "lightblue",
                "Plasma" = "royalblue4",
                "DC" = "sienna4",
                "DC2" = "sienna3",
                "pDC" = "sienna1",
                "GMP" = "orange",
                "Mast" = "brown4",
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
                "Cycling T" = "#490092",
                "CD16+ NK" = "aquamarine",
                "CD16- NK" = "aquamarine3",
                "NK" = "aquamarine4",
                "ILC" = "lightseagreen",
                "Early MK" = "gray75",
                "HSC/MPP" = "gray50")

col_order <- c("Memory B cells", "Naive B cells", 
               "Plasma cells",
               "DC", "DC2",
               "pDC", "GMP",
               "Mast cells",
               "Classical monocytes", "Non-classical monocytes",
               "Macrophages",
               "Follicular helper T cells", "MAIT cells", "Regulatory T cells",
               "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells",
               "Tem/Effector cytotoxic T cells", "Tem/Effector helper T cells",
               "Cycling T cells",
               "CD16+ NK cells", "CD16- NK cells", "NK cells", "ILC",
               "Early MK",
               "HSC/MPP")

SE$subclass <- gsub(" cells", "", SE$subclass)

SE_plot <- ggplot(SE, aes(x = npass, y = SE, color = subclass, group = subclass)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = col_scheme,
                    breaks = gsub(" cells", "", col_order),
                    name = "")+
  theme_classic() +
  theme(axis.text = element_text(color = "black"))

comp <- data.frame("comp" = c(fit$subclass$comp_amount,
                            fit2$subclass$comp_amount,
                            fit3$subclass$comp_amount,
                            fit5$subclass$comp_amount),
                   "subclass" = names(fit$subclass$se),
                   "npass" = rep(c(1, 2, 3, 4),
                                 each = length(fit$subclass$se)))

comp$subclass <- gsub(" cells", "", comp$subclass)

comp_plot <- ggplot(comp, aes(x = npass, y = comp, color = subclass, group = subclass)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = col_scheme,
                     breaks = gsub(" cells", "", col_order),
                     name = "")+
  theme_classic() +
  theme(axis.text = element_text(color = "black"))

pdf("SE_comp_typist_peac_nolegend.pdf", width = 5, height = 2.5)
ggarrange(SE_plot + theme(legend.position = "none"),
          comp_plot + theme(legend.position = "none"))
dev.off()

####CZ gene bubble plot####

cz_gene <- read.csv("CELLxGENE_gene_expression_080425.csv", skip = 9)
cz_gene_bld <- cz_gene[cz_gene$Tissue == "blood", ]
cz_gene_bld$percent_expressed <- (cz_gene_bld$Number.of.Cells.Expressing.Genes/cz_gene_bld$Cell.Count) * 100

genes <- c("HBA1", "HBB", 
           "IFIT2", "IFIT3", "OAS3", "VNN2", "AQP9", "LTF", "FFAR2", "DYSF",
           "PTGS2", "FPR2",
           "DEFA3")

setdiff(names(fit5$subclass$removed), genes)

cells <- c("granulocyte", "neutrophil", "immature neutrophil",
           "basophil", "erythroid lineage cell", "enucleated reticulocyte",
           "erythrocyte", "erythroblast", "blood cell")

cz_gene_bld_sub <- cz_gene_bld[cz_gene_bld$Gene.Symbol %in% genes &
                                 cz_gene_bld$Cell.Type %in% cells, ]

cz_gene_bld_sub$Gene.Symbol <- factor(cz_gene_bld_sub$Gene.Symbol,
                                      levels = genes)

cz_gene_bld_sub$Cell.Type <- factor(cz_gene_bld_sub$Cell.Type,
                                    levels = rev(cells))

library(scales)

pdf("cz_bubble_plot.pdf", width = 6, height = 4)
ggplot(cz_gene_bld_sub, 
       aes(x = Gene.Symbol, y = Cell.Type, fill = Expression..Scaled,
           size = percent_expressed)) +
  geom_point(pch = 21) +
  scale_fill_viridis_c(direction = -1, 
                       guide = guide_colorbar(raster = FALSE)) +
  labs(x = "", y = "", size = "Expressed in \ncells (%)",
       fill = "Expression") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust = 0.5))
dev.off()





