# load packages -----------------------------------------------------------

library(magrittr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pROC)
library(igraph)
library(RColorBrewer)
library(ggalluvial)
library(pheatmap)
source("final-code/analysis_data_functions.R", encoding = "UTF-8")

set.seed(38)

sce_ad <- Read10X_h5("data/GSE129308/GSM3704357_1-MAP2_filtered_feature_bc_matrix.h5") %>% 
  data.frame() %>% 
  xxs_filter_matrix_by_gene(gene = data_co_gene)
sce_nc <- Read10X_h5("data/GSE129308/GSM6261344_Control-1-MAP2_filtered_feature_bc_matrix.h5") %>% 
  data.frame() %>% 
  xxs_filter_matrix_by_gene(gene = data_co_gene)

set.seed(38)
train_index_ad <- sample(1:ncol(sce_ad),  ncol(sce_ad) * 0.7)
train_index_nc <- sample(1:ncol(sce_nc),  ncol(sce_nc) * 0.7)
sce_ad_test <- sce_ad[, -train_index_ad]
sce_nc_test <- sce_nc[, -train_index_nc]
sce_ad <- sce_ad[, train_index_ad]
sce_nc <- sce_nc[, train_index_nc]

# 内部拆分为训练和验证组 -------------------------------------------------------------

pair_data_sce <- xxs_get_pair_data_from_matrix(sce_ad, sce_nc, rate_cutoff = 0.6)
lasso_pair_sce <- xxs_get_pair_from_lasso(pair_data_sce, sce_ad, sce_nc)

lasso_pair_plot <- subset(lasso_pair_sce, subset = lasso_pair_sce$Freq > 90)
lasso_pair_plot <- lasso_pair_plot[order(lasso_pair_plot$Freq, decreasing = F), ]
lasso_pair_plot$. <- factor(lasso_pair_plot$., levels = lasso_pair_plot$.)
colnames(lasso_pair_plot) <- c("Gene Pair", "Freq")
ggplot(lasso_pair_plot[1:10, ]) +
  geom_col(aes(x = Freq, y = `Gene Pair`)) +
  theme_classic()

lasso_pair_sce_slc <- subset(lasso_pair_sce, subset = lasso_pair_sce$Freq > 90, select = .) %>% 
  .[, 1] %>% 
  str_split(., "_vs_", simplify = T)

sce_label_train <- c(rep(1, ncol(sce_ad)), rep(0, ncol(sce_nc))) 
sce_matrix_cbind_train <- xxs_get_adj_matrix(lasso_pair_sce_slc, sce_ad, sce_nc)

sce_cv.fit <- cv.glmnet(sce_matrix_cbind_train, sce_label_train, family = "binomial")

Lambda <- sce_cv.fit$lambda.1se
Coefficients <- coef(sce_cv.fit, s = Lambda)
coeff <- as.matrix(Coefficients)[-1]
Active.Index <- which(coeff != 0)
Active.Coefficients_sce <- coeff[Active.Index]
names(Active.Coefficients_sce) <- row.names(Coefficients)[-1][Active.Index]
Active.Coefficients_sce <- data.frame(Active.Coefficients_sce)
rank <- order(Active.Coefficients_sce$Active.Coefficients_sce, decreasing = F)
Active.Coefficients_sce$`Gene Pair` <- factor(row.names(Active.Coefficients_sce), levels = row.names(Active.Coefficients_sce)[rank])
print(Active.Coefficients_sce)

ggplot(Active.Coefficients_sce[1:10, ]) +
  geom_point(aes(x = Active.Coefficients_sce, y = `Gene Pair`)) +
  theme(axis.title.y = NA) +
  theme_classic()

pred_train <- predict(sce_cv.fit, sce_matrix_cbind_train, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
auc1 <- roc(sce_label_train, pred_train, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Training")


# sce_ad_test <- Read10X_h5("data/GSE129308/GSM3704367_6-MAP2_filtered_feature_bc_matrix.h5") %>% 
#   data.frame() %>% 
#   xxs_filter_matrix_by_gene(gene = data_co_gene)
# sce_nc_test <- Read10X_h5("data/GSE129308/GSM6261349_Control-6-MAP2_filtered_feature_bc_matrix.h5") %>% 
#   data.frame() %>% 
#   xxs_filter_matrix_by_gene(gene = data_co_gene)

sce_label_test <- c(rep(1, ncol(sce_ad_test)), rep(0, ncol(sce_nc_test))) 
sce_matrix_cbind_test <- xxs_get_adj_matrix(lasso_pair_sce_slc, sce_ad_test, sce_nc_test)

pred_test <- predict(sce_cv.fit, sce_matrix_cbind_test, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
auc2 <- roc(sce_label_test, pred_test, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Testing")

# validation --------------------------------------------------------------

xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "GSE48350", 
                        "data/GSE48350/GSE48350_series_matrix_ad.csv",
                        "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "GSE5281",
                        "data/GSE5281/GSE5281_series_matrix_ad.csv",
                        "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "GSE33000",
                        "data/GSE33000/GSE33000_series_matrix_ad.csv",
                        "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "GSE104704", 
                        "data/GSE104704/GSE104704_RNAseq_ad.csv",
                        "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "GSE159699", 
                        "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                        "data/GSE159699/GSE159699_RNAseq_count_nc.csv")
