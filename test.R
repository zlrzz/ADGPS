auc_all <- matrix(0, 100, 14) %>% 
  data.frame()
colnames(auc_all) <- c("bulk_train", "bulk_test", "b1", "b2", "b3", "b4", "b5", 
                       "sce_train", "sce_test", "s1", "s2", "s3", "s4", "s5")
for(i in 1:100) {
set.seed(i)
train_index_ad <- sample(1:ncol(AD),  ncol(AD) * 0.7)
train_index_nc <- sample(1:ncol(NC),  ncol(NC) * 0.7)
matrix_train_ad <- AD[, train_index_ad]
matrix_train_nc <- NC[, train_index_nc]
matrix_test_ad <- AD[, -train_index_ad]
matrix_test_nc <- NC[, -train_index_nc]

pair_data_from_bulk <- xxs_get_pair_data_from_matrix(matrix_train_ad, matrix_train_nc, rate_cutoff = 0.6)

lasso_pair_from_bulk <- xxs_get_pair_from_lasso(pair_data_from_bulk, matrix_train_ad, matrix_train_nc)
lasso_pair_from_sce <- xxs_get_pair_from_lasso(pair_data_from_sce, matrix_train_ad, matrix_train_nc)

lasso_pair_slc_from_bulk <- subset(lasso_pair_from_bulk, subset = lasso_pair_from_bulk$Freq > 60, select = .) %>% 
  .[, 1] %>% 
  str_split(., "_vs_", simplify = T)
lasso_pair_slc_from_sce <- subset(lasso_pair_from_sce, subset = lasso_pair_from_sce$Freq > 60, select = .) %>% 
  .[, 1] %>% 
  str_split(., "_vs_", simplify = T)

label_train <- c(rep(1, ncol(matrix_train_ad)), rep(0, ncol(matrix_train_nc))) 
label_test <- c(rep(1, ncol(matrix_test_ad)), rep(0, ncol(matrix_test_nc))) 

matrix_cbind_train_from_bulk <- xxs_get_adj_matrix(lasso_pair_slc_from_bulk, matrix_train_ad, matrix_train_nc)
matrix_cbind_test_from_bulk <- xxs_get_adj_matrix(lasso_pair_slc_from_bulk, matrix_test_ad, matrix_test_nc)
matrix_cbind_train_from_sce <- xxs_get_adj_matrix(lasso_pair_slc_from_sce, matrix_train_ad, matrix_train_nc)
matrix_cbind_test_from_sce <- xxs_get_adj_matrix(lasso_pair_slc_from_sce, matrix_test_ad, matrix_test_nc)

cv.fit_from_bulk <- cv.glmnet(matrix_cbind_train_from_bulk, label_train, family = "binomial")
cv.fit_from_sce <- cv.glmnet(matrix_cbind_train_from_sce, label_train, family = "binomial")

pred_train_from_bulk <- predict(cv.fit_from_bulk, matrix_cbind_train_from_bulk, type = "response") %>%
  as.numeric()
auc1 <- roc(label_train, pred_train_from_bulk, levels = c("0", "1"), direction = "<")
pred_test_from_bulk <- predict(cv.fit_from_bulk, matrix_cbind_test_from_bulk, type="response") %>%  #s=seLambda
  as.numeric()
auc2 <- roc(label_test, pred_test_from_bulk, levels = c("0", "1"),direction = "<")

auc_bulk_train <- round(auc1$auc, 4)
auc_bulk_test <- round(auc2$auc, 4)

pred_train_from_sce <- predict(cv.fit_from_sce, matrix_cbind_train_from_sce, type = "response") %>%
  as.numeric()
auc1 <- roc(label_train, pred_train_from_sce, levels = c("0", "1"), direction = "<")
pred_test_from_sce <- predict(cv.fit_from_sce, matrix_cbind_test_from_sce, type="response") %>%  #s=seLambda
  as.numeric()
auc2 <- roc(label_test, pred_test_from_sce, levels = c("0", "1"),direction = "<")

auc_sce_train <- round(auc1$auc, 4)
auc_sce_test <- round(auc2$auc, 4)

test_1 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_bulk, cv.fit_from_bulk, "GSE48350", 
                        "data/GSE48350/GSE48350_series_matrix_ad.csv",
                        "data/GSE48350/GSE48350_series_matrix_nc.csv")
test_2 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_bulk, cv.fit_from_bulk, "GSE5281",
                        "data/GSE5281/GSE5281_series_matrix_ad.csv",
                        "data/GSE5281/GSE5281_series_matrix_nc.csv")
test_3 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_bulk, cv.fit_from_bulk, "GSE33000",
                        "data/GSE33000/GSE33000_series_matrix_ad.csv",
                        "data/GSE33000/GSE33000_series_matrix_nc.csv")
test_4 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_bulk, cv.fit_from_bulk, "GSE104704", 
                        "data/GSE104704/GSE104704_RNAseq_ad.csv",
                        "data/GSE104704/GSE104704_RNAseq_nc.csv")
test_5 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_bulk, cv.fit_from_bulk, "GSE159699", 
                        "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                        "data/GSE159699/GSE159699_RNAseq_count_nc.csv")

test_6 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_sce, cv.fit_from_sce, "GSE48350", 
                        "data/GSE48350/GSE48350_series_matrix_ad.csv",
                        "data/GSE48350/GSE48350_series_matrix_nc.csv")
test_7 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_sce, cv.fit_from_sce, "GSE5281",
                        "data/GSE5281/GSE5281_series_matrix_ad.csv",
                        "data/GSE5281/GSE5281_series_matrix_nc.csv")
test_8 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_sce, cv.fit_from_sce, "GSE33000",
                        "data/GSE33000/GSE33000_series_matrix_ad.csv",
                        "data/GSE33000/GSE33000_series_matrix_nc.csv")
test_9 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_sce, cv.fit_from_sce, "GSE104704", 
                        "data/GSE104704/GSE104704_RNAseq_ad.csv",
                        "data/GSE104704/GSE104704_RNAseq_nc.csv")
test_10 <- xxs_auc_from_lasso_pair(lasso_pair_slc_from_sce, cv.fit_from_sce, "GSE159699", 
                        "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                        "data/GSE159699/GSE159699_RNAseq_count_nc.csv")
auc_all[i, ] <- c(auc_bulk_train, auc_bulk_test, test_1, test_2, test_3, test_4, test_5, 
                  auc_sce_train, auc_sce_test, test_6, test_7, test_8, test_9, test_10)
}

write.csv(auc_all, "results/auc_all.csv")
