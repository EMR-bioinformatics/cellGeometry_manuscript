library(rpart)
library(groupdata2)
library(rpart.plot)

load("/media/gcpeac/Rachel/Packages/cellGeometry_paper/AMP_pathotype/Combine/cellgeo_outputs.rdata")

# Test one repeat
outputdf_noNA <- outputdf[is.na(outputdf$Pathotype_orig) == FALSE, ]
outputdf_noNA$Lymphoid <- NA
outputdf_noNA$Lymphoid <- ifelse(outputdf_noNA$Pathotype_orig == "Lymphoid", "Lymphoid", "Non_lymphoid")

fit_lym <- rpart(Lymphoid ~ B_cell, data = outputdf_noNA)

outputdf_Blow <- outputdf[outputdf$B_cell > 30, ]
outputdf_Blow$Myeloid <- NA
outputdf_Blow$Myeloid <- ifelse(outputdf_Blow$Pathotype_orig == "Myeloid", "Myeloid", "Non_myeloid")
fit_mye <- rpart(Myeloid ~ Monocyte, data = outputdf_Blow, maxdepth = 1,
                 method = 'class')

par(mfrow = c(1, 2))
rpart.plot(fit_lym,
           box.palette = list(alpha("steelblue1", 0.7),
                              alpha("green3", 0.7)),
           split.font = 1,
           legend.cex = 0)
rpart.plot(fit_mye,
           box.palette = list(alpha("firebrick1", 0.7),
                              alpha("green3", 0.7)),
           split.font = 1,
           legend.cex = 0)

#9/10 folds with repeats
tree_res2step <- lapply(1:15, function(x){
  folddf <- fold(outputdf_noNA, k = 10, cat_col = "Lymphoid")
  
  fit_lym_train <- rpart(Lymphoid ~ B_cell, data = folddf[folddf$.folds %in% c(1:9), ],
                         maxdepth = 1,
                         method = 'class')
  fit_lym_test <- rpart(Lymphoid ~ B_cell, data = folddf[folddf$.folds == 10, ],
                        maxdepth = 1,
                        method = 'class')
  
  folddfv2 <- fold(outputdf_Blow, k = 10, cat_col = "Myeloid")
  
  fit_mye_train <- rpart(Myeloid ~ Monocyte, data = folddfv2[folddfv2$.folds %in% c(1:9), ],
                         maxdepth = 1,
                         method = 'class')
  fit_mye_test <- rpart(Myeloid ~ Monocyte, data = folddfv2[folddfv2$.folds == 10, ],
                        maxdepth = 1,
                        method = 'class')
  
  list(lym_folds = folddf$.folds,
       mye_folds = folddfv2$.folds,
       lym_train = fit_lym_train,
       lym_test = fit_lym_test,
       mye_train = fit_mye_train,
       mye_test = fit_mye_test)
})

names(tree_res2step) <- paste0("Rep", 1:15)

#saveRDS(tree_res2step, file = "rpart_thresholds2step.rds")

tree_res2stepdf <- data.frame()

for(i in names(tree_res2step)){
  temp <- c("B_cell_train" = tree_res2step[[i]]$lym_train[["splits"]][4],
            "B_cell_test" = tree_res2step[[i]]$lym_test[["splits"]][4],
            "Monocyte_train" = tree_res2step[[i]]$mye_train[["splits"]][4],
            "Monocyte_test" = ifelse(is.null(tree_res2step[[i]]$mye_test[["splits"]][4]),
                                     NA,
                                     tree_res2step[[i]]$mye_test[["splits"]][4]))
  
  tree_res2stepdf <- rbind(tree_res2stepdf,
                           temp)
}

tree_res2stepdf <- as.data.frame(t(tree_res2stepdf))

tree_res2stepdf$Celltype <- rep(c("B cell", "Monocyte"), each = 2)
tree_res2stepdf$fold <- rep(c("9 folds", "1 fold"), 2)

tree_res2stepdf_long <- pivot_longer(tree_res2stepdf,
                                     cols = 1:15,
                                     names_to = "repeat",
                                     values_to = "threshold")


tree_res2stepdf_long$fold <- factor(tree_res2stepdf_long$fold,
                                    levels = c("9 folds",
                                               "1 fold"))

tree_res2stepdf_long$threshold <- as.numeric(tree_res2stepdf_long$threshold)    

tree_res2stepdf_long %>% 
  group_by(fold, Celltype) %>%
  summarise(mean = mean(threshold, na.rm = TRUE),
            gmean = Gmean(threshold, na.rm = TRUE))
# fold    Celltype  mean gmean
# <fct>   <chr>    <dbl> <dbl>
#   1 9 folds B cell    25.0  24.5
# 2 9 folds Monocyte 388.  387. 
# 3 1 fold  B cell    36.2  30.5
# 4 1 fold  Monocyte 437.  419. 

tree_res2stepdf_long$Celltype[tree_res2stepdf_long$Celltype == "Monocyte"] <- "Macrophage"

ggplot(tree_res2stepdf_long[tree_res2stepdf_long$fold == "9 folds", ], 
       aes(x = Celltype, y = threshold)) +
  geom_jitter(width = 0.2, alpha = 0.8) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  scale_y_continuous(trans='log10') +
  facet_wrap(~Celltype, scale = "free") +
  labs(x = "", y = "Threshold") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank())