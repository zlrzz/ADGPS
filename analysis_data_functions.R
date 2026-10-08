# functions for analysis data 

# packages ----------------------------------------------------------------

library(glmnet)
library(RCircos)
library(magrittr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pROC)
library(igraph)
library(RColorBrewer)
library(Seurat)
library(glmGamPoi)
library(MetBrewer)
library(VennDiagram)
library(ggalluvial)
library(pheatmap)
library(ggExtra)

# xxs_filter_matrix_by_gene -----------------------------------------------

xxs_filter_matrix_by_gene <- function(matrix, gene) {
  index <- na.omit(match(gene, rownames(matrix)))
  matrix <- matrix[index, ]
  return(matrix)
}

# xxs_get_pair_data_from_matrix -------------------------------------------

xxs_get_pair_data_from_matrix <- function(matrix_A, matrix_B, p.adj_cutoff = 0.05, rate_cutoff = 0.6) { # A：疾病；B：健康
  
  co_gene <- intersect(rownames(matrix_A), rownames(matrix_B)) # 得到疾病和健康两组矩阵的的公共基因
  matrix_A <- matrix_A[match(co_gene, rownames(matrix_A)), ] # 整理矩阵，使基因一致
  matrix_B <- matrix_B[match(co_gene, rownames(matrix_B)), ]
  gene_pair <- data.frame(expand.grid(co_gene, co_gene))[, c(2, 1)] # 构建基因对 匹配后续方程  换个位置
  
  library(parallel) # 并行包
  all_cores <- detectCores() # 
  used_cores <- makeCluster(all_cores - 1) # 不讲道理直接用的剩一个 
  
  compare_pairs <- function(input_matrix) { # 从矩阵数据得到配对大小的比较总和
    num_1 <- as.numeric(parApply(used_cores, input_matrix, 1, function(x) 
    {apply(input_matrix, 1, function(y) {sum(x > y, na.rm = T)})}))
    num_2 <- as.numeric(parApply(used_cores, input_matrix, 1, function(x) 
    {apply(input_matrix, 1, function(y) {sum(x <= y, na.rm = T)})}))
    outcomes <- data.frame(num_1, num_2)
    return(outcomes)
  }
  
  num <- data.frame(compare_pairs(matrix_A), compare_pairs(matrix_B))
  pair_data <- data.frame(gene_pair, num) # 整理数据
  
  # 计算相对占比，该比率为gene-A为疾病阳性的判别比例
  pair_data$rate <- parApply(used_cores, pair_data, 1, function(x) {
    num_4 <- as.numeric(x[c(3, 4, 5, 6)])
    rate <- (num_4[1] + num_4[4]) / sum(num_4)
    return(rate)
  })
  
  unique_pair <- function(pair_data) { # 去除重复基因对（选效果好的）
    rate <- abs(pair_data$rate - 0.5)
    rate_matrix <- matrix(rate, nrow = length(co_gene), byrow = T)
    n <- nrow(rate_matrix)
    index_slc <- c()
    for(i in 1:n) {
      for(j in 1:n) {
        if(i >= j) {
          next
        }
        index <- ifelse(rate_matrix[i, j] >= rate_matrix[j, i], i*n - n + j, j*n - n + i)
        index_slc <- c(index_slc, index)
      }
    }
    output_data <- pair_data[index_slc, ]
    return(output_data)
  }
  
  pair_data <- unique_pair(pair_data)
  
  # 将基因对更改为前者为疾病阳性基因，后者为疾病阴性基因。并按照比率进行排序
  index <- which(pair_data$rate < 0.5)
  pair_data[index, ] <- pair_data[index, c(2, 1, 4, 3, 6, 5, 7)]
  pair_data[index, 7] <- 1 - pair_data[index, 7]
  pair_data <- pair_data[order(pair_data$rate, decreasing = T), ]
  
  p.value <- parApply(used_cores, pair_data, 1, function(x) {
    num_5 <- as.numeric(x[c(3, 4, 5, 6)])
    fisher.test(matrix(num_5, 2, 2))$p.value
  })
  p.adj <- p.adjust(p.value, method = "BH")
  
  pair_data <- data.frame(pair_data, p.adj)
  colnames(pair_data) <- c("Gene-A", "Gene-B", 
                           "A>B in D", "A<B in D", 
                           "A>B in N", "A<B in N", 
                           "rate", "p.adj")
  pair_data$gene_pair <- paste0(pair_data[, 1], "_vs_", pair_data[, 2])
  pair_data <- subset(pair_data, pair_data$p.adj < p.adj_cutoff & pair_data$rate > rate_cutoff) # 去除掉p小于0.05
  stopCluster(used_cores)
  return(pair_data) # 返回数据
}

# xxs_step_forward --------------------------------------------------------

xxs_step_forward <- function(pair_data, matrix_A, matrix_B, step_length = 0) {
  
  library(pROC)
  label <- c(rep(1, ncol(matrix_A)), rep(0, ncol(matrix_B))) # 样本信息
  co_gene <- intersect(rownames(matrix_A), rownames(matrix_B))
  all_matrix <- cbind(matrix_A[match(co_gene, rownames(matrix_A)), ],
                      matrix_B[match(co_gene, rownames(matrix_B)), ])
  steps <- c(1:nrow(pair_data))
  scores <- 0
  slc_pair <- data.frame()
  max_auc <- 0
  
  for (step in steps) {
    # 每步的打分
    score_A <- as.numeric(all_matrix[match(pair_data[step, 1], co_gene), ])
    score_B <- as.numeric(all_matrix[match(pair_data[step, 2], co_gene), ])
    score <- score_A - score_B
    score[score > 0] <- 1
    score[score < 0] <- 0
    # 该步的roc
    pred_socres <- (score + scores) / (nrow(slc_pair) + 1)
    roc_auc <- round(roc(label, pred_socres, levels = c("0", "1"), direction = "<")$auc, 4)
    # 提高入选
    if (roc_auc > (max_auc + step_length)) {
      max_auc <- roc_auc
      # print(max_auc)
      pair_data[step, 10] <- max_auc
      slc_pair <- rbind(slc_pair, pair_data[step, ])
      scores <- score + scores
    }
  }
  
  return(slc_pair)
  
}

# xxs_get_pair_from_lasso -------------------------------------------------

xxs_get_pair_from_lasso <- function(pair_data, matrix_a, matrix_b, lasso_times = 100) {
  library(glmnet)
  library(parallel) # 并行包
  all_cores <- detectCores() # 
  used_cores <- makeCluster(all_cores - 1) # 不讲道理直接用的剩一个
  # 为lasso调整数据
  label <- c(rep(1, ncol(matrix_a)), rep(0, ncol(matrix_b))) # 正常样本标记为0,疾病样本标记为1
  matrix_cbind <- cbind(matrix_a, matrix_b) # 表达矩阵整合
  clusterExport(used_cores, c("matrix_cbind", "pair_data"),envir = environment())
  matrix_cbind <- parApply(used_cores, pair_data, 1, function(x) { # 根据配对信息转化为0/1矩阵
    library(magrittr)
    gene_a_exp <- x[1] %>% 
      match(., rownames(matrix_cbind)) %>% 
      matrix_cbind[., ] %>% 
      as.numeric()
    gene_b_exp <- x[2] %>% 
      match(., rownames(matrix_cbind)) %>% 
      matrix_cbind[., ] %>% 
      as.numeric()
    num <- gene_a_exp - gene_b_exp
    num[num > 0] <- 1
    num[num < 0] <- 0
    num[is.na(num)] <- 0
    num <- as.numeric(num)
    return(num)
  })
  stopCluster(used_cores) # 关闭集群
  
  colnames(matrix_cbind) <- pair_data$gene_pair
  
  for (p in 1:100) {
    
    set.seed(p)
    train_index <- sample(1:nrow(matrix_cbind),  nrow(matrix_cbind) * 0.7)
    matrix_train <- matrix_cbind[train_index, ]
    label_train <- label[train_index]
    cv.fit <- cv.glmnet(matrix_train, label_train, family = "binomial")
    coef_of_pair <- cv.fit$lambda.1se %>%  # 1.得到此次一个标准差内的最大lambda
      coef(cv.fit, s = .) %>% # 2.模型系数
      as.matrix() %>% # 3.格式转换
      .[-1, ] # 4. 整理数据
    lasso_pair_once <- which(coef_of_pair != 0) %>% # 不等于0的系数
      names() # 提取名字
    if(p == 1) {
      lasso_pair <- lasso_pair_once
      next
    }
    lasso_pair <- c(lasso_pair, lasso_pair_once)
  }
  
  lasso_pair %>% 
    table() %>% 
    as.data.frame() %>% 
    return()
}

# xxs_auc_from_lasso_pair -------------------------------------------------

xxs_auc_from_lasso_pair <- function(lasso_pair_slc, cv.fit, save_dir, title, dir_a, dir_b) {
  
  slc_gene <- unique(as.character(lasso_pair_slc[, c(1, 2)]))
  matrix_a <- read.csv(dir_a, row.names = 1) %>% xxs_filter_matrix_by_gene(., slc_gene)
  matrix_b <- read.csv(dir_b, row.names = 1) %>% xxs_filter_matrix_by_gene(., slc_gene)
  
  label <- c(rep(1, ncol(matrix_a)), rep(0, ncol(matrix_b))) # 正常样本标记为0,癌症样本标记为1
  matrix_cbind <- cbind(matrix_a, matrix_b)  # 表达矩阵整合
  matrix_cbind <- apply(lasso_pair_slc, 1, function(x) { # 根据配对信息转化为0/1矩阵
    gene_a_exp <- x[1] %>% 
      match(., rownames(matrix_cbind)) %>% 
      matrix_cbind[., ] %>% 
      as.numeric()
    gene_b_exp <- x[2] %>% 
      match(., rownames(matrix_cbind)) %>% 
      matrix_cbind[., ] %>% 
      as.numeric()
    num <- gene_a_exp - gene_b_exp
    num[num > 0] <- 1
    num[num < 0] <- 0
    num <- as.numeric(num)
    return(num)
  })
  
  pred_test <- predict(cv.fit, matrix_cbind, type = "response") %>%  #s=seLambda
    as.numeric()
  png(paste0(save_dir, title, " test.png"), width = 900, height = 900, units = "px", res = 300)
  auc <- roc(label, pred_test, main = title,
             levels = c("0", "1"), plot = T, print.auc = T, direction = "<", legacy.axes = T)
  dev.off()
  print(round(auc$auc, 4))
  
}

# xxs_get_adj_matrix ------------------------------------------------------

xxs_get_adj_matrix <- function(pair_slc, matrix_a, matrix_b) {
  matrix_cbind_ori <- cbind(matrix_a, matrix_b) # 表达矩阵整合
  matrix_cbind <- apply(pair_slc, 1, function(x) { # 根据配对信息转化为0/1矩阵
    gene_a_exp <- x[1] %>% 
      match(., rownames(matrix_cbind_ori)) %>% 
      matrix_cbind_ori[., ] %>% 
      as.numeric()
    gene_b_exp <- x[2] %>% 
      match(., rownames(matrix_cbind_ori)) %>% 
      matrix_cbind_ori[., ] %>% 
      as.numeric()
    num <- gene_a_exp - gene_b_exp
    num[num > 0] <- 1
    num[num < 0] <- 0
    num[is.na(num)] <- 0
    num <- as.numeric(num)
    return(num)
  })
  colnames(matrix_cbind) <- paste0(pair_slc[, 1], "_vs_", pair_slc[, 2])
  return(matrix_cbind)
}

# xxs_get_adj_matrix_1 ------------------------------------------------------

xxs_get_adj_matrix_1 <- function(pair_slc, matrix) {
  matrix <- apply(pair_slc, 1, function(x) { # 根据配对信息转化为0/1矩阵
    gene_a_exp <- x[1] %>% 
      match(., rownames(matrix)) %>% 
      matrix[., ] %>% 
      as.numeric()
    gene_b_exp <- x[2] %>% 
      match(., rownames(matrix)) %>% 
      matrix[., ] %>% 
      as.numeric()
    num <- gene_a_exp - gene_b_exp
    num[num > 0] <- 1
    num[num < 0] <- -1
    num[is.na(num)] <- 0
    num <- as.numeric(num)
    return(num)
  })
  colnames(matrix) <- paste0(pair_slc[, 1], "_vs_", pair_slc[, 2])
  return(matrix)
}

# xxs_exp_col_plot --------------------------------------------------------

xxs_exp_col_plot <- function(lasso_pair_slc, dir, title, dir_ad, dir_nc) {
  EXP_ad <- read.csv(dir_ad, row.names = 1) %>% xxs_get_adj_matrix_1(lasso_pair_slc, .)
  EXP_nc <- read.csv(dir_nc, row.names = 1) %>% xxs_get_adj_matrix_1(lasso_pair_slc, .)
  
  EXP_ad <- apply(EXP_ad, 2, function(x) sum(x == 1)/nrow(EXP_ad) * 100) %>% data.frame()
  EXP_ad$Gene_pair <- rownames(EXP_ad)
  EXP_nc <- apply(EXP_nc, 2, function(x) sum(x == 1)/nrow(EXP_nc) * 100) %>% data.frame()
  EXP_nc$Gene_pair <- rownames(EXP_nc)
  
  EXP <- rbind(EXP_ad, EXP_nc)
  n <- nrow(lasso_pair_slc)
  EXP$Group <- c(rep("AD", n), rep("NC", n))
  EXP$Gene_pair <- factor(EXP$Gene_pair, levels = rownames(EXP_ad))
  ggplot(EXP) +
    geom_col(aes(y = Gene_pair, x = ., fill = Group), width = 0.6, position = position_dodge2()) +
    theme_classic() 
  # +
  # xlab("Percent (%)") +
  # ggtitle(title)
  ggsave(paste0(dir, title, " barplot.png"), width = 4, height = 6, dpi = 300)
}

# xxs_boxplot -------------------------------------------------------------

xxs_boxplot <- function(lasso_pair_slc, dir, title, dir_ad, dir_nc) {
  gene_plot <- c(lasso_pair_slc) %>% 
    unique()
  EXP_ad <- read.csv(dir_ad, row.names = 1) %>% 
    xxs_filter_matrix_by_gene(., gene_plot)
  EXP_nc <- read.csv(dir_nc, row.names = 1) %>% 
    xxs_filter_matrix_by_gene(., gene_plot)
  exp <- cbind(EXP_ad, EXP_nc)
  exp[is.na(exp)] <- 0 
  
  if(sum(as.numeric(unlist(exp)) > 100) > 0) {
    exp <- log2(exp + 1)
  } 
  
  if(sum(as.numeric(unlist(exp)) < 0) == 0) {
    exp_mean <- apply(exp, 2, mean)
    sd_mean <- apply(exp, 2, sd)
    exp <- (exp - exp_mean) / sd_mean
  }
  
  for(i in 1:nrow(lasso_pair_slc)) {
    genea <- lasso_pair_slc[i, 1]
    geneb <- lasso_pair_slc[i, 2]
    expa <- which(rownames(EXP_ad) == genea) %>% 
      exp[., ] %>% 
      as.numeric() %>% 
      data.frame()
    expa$Gene <- genea
    expa$Group <- c(rep("AD", ncol(EXP_ad)), rep("NC", ncol(EXP_nc)))
    expb <- which(rownames(EXP_ad) == geneb) %>% 
      exp[., ] %>% 
      as.numeric() %>% 
      data.frame()
    expb$Gene <- geneb
    expb$Group <- c(rep("AD", ncol(EXP_ad)), rep("NC", ncol(EXP_nc)))
    exp_plot <- rbind(expa, expb)
    exp_plot$Gene <- factor(exp_plot$Gene, levels = c(genea, geneb))
    ggplot(exp_plot) +
      geom_boxplot(aes(x = Group, y = ., fill = Gene)) +
      theme_bw() +
      scale_fill_manual(values = c("goldenrod1", "dodgerblue1")) +
      ylab("Exp") +
      xlab(NULL) +
      theme(legend.position = "top", )
    ggsave(paste0(dir, title, genea, geneb, " boxplot.png"), width = 3, height = 2, dpi = 300)
  }
  
}

# xxs_waterfall -----------------------------------------------------------

xxs_waterfall_6 <- function(lasso_pair_slc, title, dir, dir_ad, dir_nc) {
  EXP_ad <- read.csv(dir_ad, row.names = 1) %>% 
    xxs_get_adj_matrix_1(lasso_pair_slc, .) %>% 
    data.frame()
  
  EXP_nc <- read.csv(dir_nc, row.names = 1) %>% 
    xxs_get_adj_matrix_1(lasso_pair_slc, .) %>% 
    data.frame()
  
  EXP_ad <- EXP_ad[order(EXP_ad[, 1],
                         EXP_ad[, 2],
                         EXP_ad[, 3],
                         EXP_ad[, 4],
                         EXP_ad[, 5],
                         EXP_ad[, 6],
                         # EXP_ad[, index[7]],
                         decreasing = T), ]
  
  EXP_nc <- EXP_nc[order(EXP_nc[, 1],
                         EXP_nc[, 2],
                         EXP_nc[, 3],
                         EXP_nc[, 4],
                         EXP_nc[, 5],
                         EXP_nc[, 6],
                         # EXP_nc[, index[7]],
                         decreasing = T), ]
  EXP_GSE48350 <- rbind(EXP_ad, EXP_nc) %>% 
    t() %>% 
    data.frame()
  
  annotation_col <- c(rep("AD", nrow(EXP_ad)),
                      rep("NC", nrow(EXP_nc))) %>% 
    data.frame()
  rownames(annotation_col) <- colnames(EXP_GSE48350)
  colnames(annotation_col) <- "Group"
  annotation_col$Group <- factor(annotation_col$Group, levels = c("AD", "NC"))
  anno_col <- list(Group = c(AD = '#F8766D', NC = '#00BFC4')) 
  png(paste0(dir, title, " heatmap.png"), width = 1250, height = 750, res = 300)
  pheatmap(EXP_GSE48350, 
           gaps_col = nrow(EXP_ad),
           annotation_col = annotation_col, 
           annotation_colors = anno_col, 
           cluster_rows = F, 
           cluster_cols = F, 
           color = c("#E0FFFF", "#FFA54F"), 
           border_color = NA, 
           show_colnames = F, 
           show_rownames = T, 
           legend = F, 
           main = title)
  dev.off()
}

xxs_waterfall_7 <- function(lasso_pair_slc, title, dir, dir_ad, dir_nc) {
  EXP_ad <- read.csv(dir_ad, row.names = 1) %>% 
    xxs_get_adj_matrix_1(lasso_pair_slc, .) %>% 
    data.frame()
  
  EXP_nc <- read.csv(dir_nc, row.names = 1) %>% 
    xxs_get_adj_matrix_1(lasso_pair_slc, .) %>% 
    data.frame()
  
  EXP_ad <- EXP_ad[order(EXP_ad[, 1],
                         EXP_ad[, 2],
                         EXP_ad[, 3],
                         EXP_ad[, 4],
                         EXP_ad[, 5],
                         EXP_ad[, 6],
                         EXP_ad[, 7],
                         decreasing = T), ]
  
  EXP_nc <- EXP_nc[order(EXP_nc[, 1],
                         EXP_nc[, 2],
                         EXP_nc[, 3],
                         EXP_nc[, 4],
                         EXP_nc[, 5],
                         EXP_nc[, 6],
                         EXP_nc[, 7],
                         decreasing = T), ]
  EXP_GSE48350 <- rbind(EXP_ad, EXP_nc) %>% 
    t() %>% 
    data.frame()
  
  annotation_col <- c(rep("AD", nrow(EXP_ad)),
                      rep("NC", nrow(EXP_nc))) %>% 
    data.frame()
  rownames(annotation_col) <- colnames(EXP_GSE48350)
  colnames(annotation_col) <- "Group"
  annotation_col$Group <- factor(annotation_col$Group, levels = c("AD", "NC"))
  anno_col <- list(Group = c(AD = '#F8766D', NC = '#00BFC4')) 
  png(paste0(dir, title, " heatmap.png"), width = 1250, height = 750, res = 300)
  pheatmap(EXP_GSE48350, 
           gaps_col = nrow(EXP_ad),
           annotation_col = annotation_col, 
           annotation_colors = anno_col, 
           cluster_rows = F, 
           cluster_cols = F, 
           color = c("#E0FFFF", "#FFA54F"), 
           border_color = NA, 
           show_colnames = F, 
           show_rownames = T, 
           legend = F, 
           main = title)
  dev.off()
}
# End ---------------------------------------------------------------------

print("functions for analysis_data load completed")
