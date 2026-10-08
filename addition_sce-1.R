# 加载包 ---------------------------------------------------------------------

source("code/analysis_data_functions.R", encoding = "UTF-8")

# integration -------------------------------------------------------------

sce_ad_seurat <- CreateSeuratObject(Read10X_h5("data/GSE129308/GSM3704357_1-MAP2_filtered_feature_bc_matrix.h5"), project = "AD")
sce_nc_seurat <- CreateSeuratObject(Read10X_h5("data/GSE129308/GSM6261344_Control-1-MAP2_filtered_feature_bc_matrix.h5"), project = "NC")
# 把你想联合分析的放一起
sce_list <- list(sce_ad_seurat, sce_nc_seurat)

# 后面都是流程
# 这四步是常用的合并的流程
# 标准化
sce_list <- lapply(X = sce_list, FUN = SCTransform)
# 找合并的基因,用前三千个高变基因
IntegrationFeatures <- SelectIntegrationFeatures(object.list = sce_list, nfeatures = 3000)
# 准备合并
sce_list <- PrepSCTIntegration(object.list = sce_list, anchor.features = IntegrationFeatures)
# 找锚定点
IntegrationAnchors <- FindIntegrationAnchors(object.list = sce_list,
                                             normalization.method = "SCT",
                                             anchor.features = IntegrationFeatures)
# 合并
combined_sce <- IntegrateData(anchorset = IntegrationAnchors, normalization.method = "SCT")
# 设置成默认的矩阵
DefaultAssay(combined_sce) <- 'integrated'

combined_sce <- RunPCA(combined_sce, verbose = FALSE)
combined_sce <- RunUMAP(combined_sce, reduction = "pca", dims = 1:30)
combined_sce <- FindNeighbors(combined_sce, reduction = "pca", dims = 1:30)
combined_sce <- FindClusters(combined_sce, resolution = 1)

# 单细胞出图以及注释 ---------------------------------------------------------------

png("results/sc/VlnPlot.orig.ident.png", width = 1200, height = 1200, res = 300)
VlnPlot(combined_sce, features = c("nCount_RNA", "nFeature_RNA"), group.by = "orig.ident", assay = "RNA", pt.size = 0)
dev.off()

Idents(combined_sce) <- combined_sce$integrated_snn_res.1
png("results/sc/DimPlot.orig.ident.png", width = 2700, height = 1200, res = 300)
DimPlot(combined_sce, reduction = "umap", split.by = "orig.ident")
dev.off()


Astrocytes <- c("GFAP", "AQP4", "SLC1A2")
Excitatory <- c("CAMK2A", "LDB2", "SLC17A7")
Inhibitory <- c("GAD1", "GAD2", "SLC32A1")
Oligodendrocytes <- c("MBP", "MOG", "PLP1")

gene_anno_list <- list(Excitatory = Excitatory,
                       Inhibitory = Inhibitory,
                       Astrocytes = Astrocytes,
                       Oligodendrocytes = Oligodendrocytes)

Idents(combined_sce) = combined_sce$integrated_snn_res.1
png("results/sc/DotPlot_gene_anno_list.png", width = 4500, height = 1500, res = 300)
DotPlot(combined_sce, features = gene_anno_list, assay = "RNA")
dev.off()


Astro <- 21
Oligo <- 20
Inhib <- c(0, 6, 10, 14, 16, 19)
Excit <- setdiff(c(1:19), Inhib)

combined_sce$seurat_clusters <- factor(combined_sce$seurat_clusters, levels = c(Astro, Oligo, Excit, Inhib))
combined_sce$celltype = ifelse(combined_sce$seurat_clusters %in% Astro, "Astrocytes",
                               ifelse(combined_sce$seurat_clusters %in% Oligo, "Oligodendrocytes",
                                      ifelse(combined_sce$seurat_clusters %in% Excit, "Excitatory",
                                             "Inhibitory")))
Idents(combined_sce) = combined_sce$celltype
combined_sce$celltype <- factor(combined_sce$celltype, levels = c("Excitatory", "Inhibitory", "Astrocytes", "Oligodendrocytes"))


png("results/sc/VlnPlot.gene_anno.png", width = 4400, height = 2200, res = 300)
VlnPlot(combined_sce, features = unlist(gene_anno_list), group.by = "celltype", assay = "RNA", pt.size = 0, ncol = 6)
dev.off()

Idents(combined_sce) = combined_sce$celltype

png("results/sc/DimPlot.celltype.png", width = 1800, height = 1200, res = 300)
DimPlot(combined_sce, reduction = "umap", label = T)
dev.off()


# load packages -----------------------------------------------------------

sce_ad <- Read10X_h5("data/GSE129308/GSM3704357_1-MAP2_filtered_feature_bc_matrix.h5") %>%
  data.frame() %>%
  xxs_filter_matrix_by_gene(gene = data_co_gene)
sce_nc <- Read10X_h5("data/GSE129308/GSM6261344_Control-1-MAP2_filtered_feature_bc_matrix.h5") %>%
  data.frame() %>%
  xxs_filter_matrix_by_gene(gene = data_co_gene)

# 内部拆分为训练和验证组 -------------------------------------------------------------

set.seed(38)
train_index_ad <- sample(1:ncol(sce_ad),  ncol(sce_ad) * 0.7)
train_index_nc <- sample(1:ncol(sce_nc),  ncol(sce_nc) * 0.7)
sce_ad_test <- sce_ad[, -train_index_ad]
sce_nc_test <- sce_nc[, -train_index_nc]
sce_ad_train <- sce_ad[, train_index_ad]
sce_nc_train <- sce_nc[, train_index_nc]

# sc模型训练 ------------------------------------------------------------------

pair_data_sce <- xxs_get_pair_data_from_matrix(sce_ad_train, sce_nc_train, rate_cutoff = 0.6)
lasso_pair_sce <- xxs_get_pair_from_lasso(pair_data_sce, sce_ad_train, sce_nc_train)
write.csv(lasso_pair_sce, "results/sc/lasso_pair_sce.csv")
lasso_pair_sce <- read.csv("results/sc/lasso_pair_sce.csv", row.names = 1)
lasso_pair_plot_sce <- subset(lasso_pair_sce, subset = lasso_pair_sce$Freq > 0)
# lasso_pair_plot_sce <- lasso_pair_plot_sce[order(lasso_pair_plot_sce$Freq, decreasing = F), ]
lasso_pair_plot_sce$. <- factor(lasso_pair_plot_sce$., levels = lasso_pair_plot_sce$.)
colnames(lasso_pair_plot_sce) <- c("Gene Pair", "Frequency")
{
  ggplot(lasso_pair_plot_sce) +
    geom_col(aes(y = -Frequency, x = `Gene Pair`, fill = Frequency)) +
    # scale_fill_gradient(low = "black", high = "#FFA54F") +
    geom_hline(yintercept = -90, col = "red", linetype = 2) +
    # guides(fill = "none") +
    coord_polar() +
    ylim(-300, NA) +
    theme_classic() +
    theme(axis.line = element_line(colour = "white")) + # 坐标轴白色
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
    xlab(NULL) + # 去掉坐标名字
    ylab(NULL) # 去掉坐标名字
  
  ggsave("results/sc/gene_pair_freq.png", device = "png", width = 8, height = 8, dpi = 300)
}


lasso_pair_sce_slc <- subset(lasso_pair_sce, subset = lasso_pair_sce$Freq > 90, select = .) %>% 
  .[, 1] %>% 
  str_split(., "_vs_", simplify = T)

sce_label_train <- c(rep(1, ncol(sce_ad_train)), rep(0, ncol(sce_nc_train))) 
sce_matrix_cbind_train <- xxs_get_adj_matrix(lasso_pair_sce_slc, sce_ad_train, sce_nc_train)

sce_cv.fit <- cv.glmnet(sce_matrix_cbind_train, sce_label_train, family = "binomial")

Lambda <- sce_cv.fit$lambda.1se
Coefficients <- coef(sce_cv.fit, s = Lambda)
coeff <- as.matrix(Coefficients)[-1]
Active.Index <- which(coeff != 0)
Active.Coefficients_sce <- coeff[Active.Index]
names(Active.Coefficients_sce) <- row.names(Coefficients)[-1][Active.Index]
Active.Coefficients_sce <- data.frame(Active.Coefficients_sce)
rank <- order(Active.Coefficients_sce$Active.Coefficients_sce, decreasing = F)
Active.Coefficients_sce$Gene_Pair <- factor(row.names(Active.Coefficients_sce), levels = row.names(Active.Coefficients_sce)[rank])
write.csv(Active.Coefficients_sce, "results/sc/Active.Coefficients_sce.csv")
Active.Coefficients_sce <- read.csv("results/sc/Active.Coefficients_sce.csv", row.names = 1)
print(Active.Coefficients_sce)
colnames(Active.Coefficients_sce)[1] <- "Coefficient"
{
  ggplot(Active.Coefficients_sce) +
    geom_col(aes(y = -abs(Coefficient), x = Gene_Pair, fill = Coefficient)) +
    coord_polar() +
    # guides(fill = "none") +
    scale_fill_gradient(low = "blue", high = "red") +
    # xlab("Active.Coefficients") +
    # theme_classic() +
    # theme(axis.line = element_line(colour = "white")) + # 坐标轴白色
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
    xlab(NULL) + # 去掉坐标名字
    ylab(NULL)  + # 去掉坐标名字 
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    ) 
  ggsave("results/sc/slc_pair_coeff.png", device = "png", width = 8, height = 8, dpi = 300)
}

pred_train <- predict(sce_cv.fit, sce_matrix_cbind_train, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
{
  png("results/sc/sc pair data train auc.png", width = 900, height = 900, units = "px", res = 300)
  auc1 <- roc(sce_label_train, pred_train, levels = c("0", "1"), legacy.axes = T, plot = T, print.auc = T, direction = "<", main = "Training")
  dev.off()
  }


sce_label_test <- c(rep(1, ncol(sce_ad_test)), rep(0, ncol(sce_nc_test))) 
sce_matrix_cbind_test <- xxs_get_adj_matrix(lasso_pair_sce_slc, sce_ad_test, sce_nc_test)

pred_test <- predict(sce_cv.fit, sce_matrix_cbind_test, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
{
  png("results/sc/sc pair data test auc.png", width = 900, height = 900, units = "px", res = 300)
  auc2 <- roc(sce_label_test, pred_test, levels = c("0", "1"), legacy.axes = T, plot = T, print.auc = T, direction = "<", main = "Testing")
  dev.off()
  }

Active.Coefficients_sce_1 <- Active.Coefficients_sce[order(abs(Active.Coefficients_sce$Coefficient), decreasing = T), ]

# validation --------------------------------------------------------------

xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "results/sc/", "GSE48350", 
                        "data/GSE48350/GSE48350_series_matrix_ad.csv",
                        "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "results/sc/", "GSE5281",
                        "data/GSE5281/GSE5281_series_matrix_ad.csv",
                        "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "results/sc/", "GSE33000",
                        "data/GSE33000/GSE33000_series_matrix_ad.csv",
                        "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "results/sc/", "GSE104704", 
                        "data/GSE104704/GSE104704_RNAseq_ad.csv",
                        "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_sce_slc, sce_cv.fit, "results/sc/", "GSE159699", 
                        "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                        "data/GSE159699/GSE159699_RNAseq_count_nc.csv")


# transfer learning -------------------------------------------------------

pair_data_from_sce <- Active.Coefficients_sce %>% 
  .[, 2] %>% 
  str_split(., "_vs_", simplify = T) %>% 
  data.frame()
pair_data_from_sce$gene_pair <- Active.Coefficients_sce$Gene_Pair

EXP_GSE5281_ad <- read.csv("data/GSE5281/GSE5281_series_matrix_ad.csv", row.names = 1) %>% xxs_filter_matrix_by_gene(., data_co_gene)
EXP_GSE5281_nc <- read.csv("data/GSE5281/GSE5281_series_matrix_nc.csv", row.names = 1) %>% xxs_filter_matrix_by_gene(., data_co_gene)

EXP_GSE48350_ad <- read.csv("data/GSE48350/GSE48350_series_matrix_ad.csv", row.names = 1) %>% xxs_filter_matrix_by_gene(., data_co_gene)
EXP_GSE48350_nc <- read.csv("data/GSE48350/GSE48350_series_matrix_nc.csv", row.names = 1) %>% xxs_filter_matrix_by_gene(., data_co_gene)

AD <- cbind(EXP_GSE5281_ad, EXP_GSE48350_ad)
NC <- cbind(EXP_GSE5281_nc, EXP_GSE48350_nc)

set.seed(38)
train_index_ad <- sample(1:ncol(AD),  ncol(AD) * 0.7)
train_index_nc <- sample(1:ncol(NC),  ncol(NC) * 0.7)
matrix_train_ad <- AD[, train_index_ad]
matrix_train_nc <- NC[, train_index_nc]
matrix_test_ad <- AD[, -train_index_ad]
matrix_test_nc <- NC[, -train_index_nc]

lasso_pair_bulk_sce_to_bulk <- xxs_get_pair_from_lasso(pair_data_from_sce, matrix_train_ad, matrix_train_nc)

# lasso_pair_plot_sce_to_bulk <- subset(lasso_pair_bulk_sce_to_bulk, subset = lasso_pair_bulk_sce_to_bulk$Freq > 30)
# lasso_pair_plot_sce_to_bulk <- lasso_pair_plot_sce_to_bulk[order(lasso_pair_plot_sce_to_bulk$Freq, decreasing = F), ]
# lasso_pair_plot_sce_to_bulk$. <- factor(lasso_pair_plot_sce_to_bulk$., levels = lasso_pair_plot_sce_to_bulk$.)
colnames(lasso_pair_bulk_sce_to_bulk) <- c("Gene Pair", "Frequency")

{
  ggplot(lasso_pair_bulk_sce_to_bulk) +
    geom_col(aes(y = -Frequency, x = `Gene Pair`, fill = Frequency)) +
    # scale_fill_gradient(low = "black", high = "#FFA54F") +
    geom_hline(yintercept = -60, col = "red", linetype = 2) +
    # guides(fill = "none") +
    coord_polar() +
    ylim(-300, NA) +
    theme_classic() +
    theme(axis.line = element_line(colour = "white")) + # 坐标轴白色
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
    xlab(NULL) + # 去掉坐标名字
    ylab(NULL) # 去掉坐标名字
  
  ggsave("results/sc_bulk/gene_pair_freq.png", device = "png", width = 8, height = 8, dpi = 300)
}

lasso_pair_slc_sce_to_bulk <- subset(lasso_pair_bulk_sce_to_bulk, subset = lasso_pair_bulk_sce_to_bulk$Freq > 60, select = `Gene Pair`) %>% 
  .[, 1] %>% 
  str_split(., "_vs_", simplify = T)

label_train <- c(rep(1, ncol(matrix_train_ad)), rep(0, ncol(matrix_train_nc))) 
matrix_cbind_train <- xxs_get_adj_matrix(lasso_pair_slc_sce_to_bulk, matrix_train_ad, matrix_train_nc)

label_test <- c(rep(1, ncol(matrix_test_ad)), rep(0, ncol(matrix_test_nc))) 
matrix_cbind_test <- xxs_get_adj_matrix(lasso_pair_slc_sce_to_bulk, matrix_test_ad, matrix_test_nc)

cv.fit <- cv.glmnet(matrix_cbind_train, label_train, family = "binomial")

Lambda <- cv.fit$lambda.1se
Coefficients <- coef(cv.fit, s = Lambda)
coeff <- as.matrix(Coefficients)[-1]
Active.Index <- which(coeff != 0)
Active.Coefficients <- coeff[Active.Index]
names(Active.Coefficients) <- row.names(Coefficients)[-1][Active.Index]
Active.Coefficients <- data.frame(Active.Coefficients)
Active.Coefficients$`Gene Pair` <- row.names(Active.Coefficients)

index <- which(Active.Coefficients$Active.Coefficients < 0)
Active.Coefficients$Active.Coefficients[index] <- -Active.Coefficients$Active.Coefficients[index]
gene_adj <- str_split(Active.Coefficients$`Gene Pair`[index], "_vs_", simplify = T)
Active.Coefficients$`Gene Pair`[index] <- paste0(gene_adj[, 2], "_vs_", gene_adj[, 1])

rank <- order(Active.Coefficients$Active.Coefficients, decreasing = F)
Active.Coefficients$`Gene Pair` <- factor(Active.Coefficients$`Gene Pair`, levels = Active.Coefficients$`Gene Pair`[rank])
Active.Coefficients$Active.Coefficients <- as.numeric(Active.Coefficients$Active.Coefficients)
print(Active.Coefficients)
write.table(Active.Coefficients, 
            "results/sc_bulk/Active.Coefficients_train_from_bulk.txt",
            row.names = F)

{ 
  ggplot(Active.Coefficients) +
    geom_col(aes(y = `Gene Pair`, x = Active.Coefficients, fill = Active.Coefficients), width = 0.5) +
    # guides(fill = guide_legend(title = 'Coefficients')) +
    # coord_polar() +
    scale_fill_gradient(low = "lightpink", high= "red3") +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      # panel.grid.major = element_blank(), # get rid of major grid
      # panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )
  ggsave("results/sc_bulk/slc_pair_coeff.png", device = "png", width = 8, height = 8, dpi = 300)
  # ggsave("results/sc_bulk/slc_pair_coeff_1.png", device = "png", width = 5, height = 3, dpi = 300)
}

pred_train <- predict(cv.fit, matrix_cbind_train, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
{
  png("results/sc_bulk/bulk pair data train auc.png", width = 900, height = 900, units = "px", res = 300)
  auc1 <- roc(label_train, pred_train, legacy.axes = T, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Training")
  dev.off()
  }# plot(auc1, main = round(auc1$auc, 4))


pred_test <- predict(cv.fit, matrix_cbind_test, type="response", s = Lambda) %>%  #s=seLambda
  as.numeric()
{
  png("results/sc_bulk/bulk pair data test auc.png", width = 900, height = 900, units = "px", res = 300)
  auc2 <- roc(label_test, pred_test, legacy.axes = T, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Testing")
  dev.off()
  }# plot(auc2, main = round(auc2$auc, 4))

ad_gene <- rownames(Active.Coefficients) %>%
  str_split(., "_vs_", simplify = T) %>%
  c() %>% 
  unique()
write.table(ad_gene, 
            "results/sc_bulk/gene_train_from_bulk.txt",
            row.names = F)
ad_gene <- read.table("results/sc_bulk/gene_train_from_bulk.txt", header = T) %>% 
  .[, 1]
# validation --------------------------------------------------------------

xxs_auc_from_lasso_pair(lasso_pair_slc_sce_to_bulk, cv.fit,
                        "results/sc_bulk/", "GSE48350", 
                        "data/GSE48350/GSE48350_series_matrix_ad.csv",
                        "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc_sce_to_bulk, cv.fit,
                        "results/sc_bulk/", "GSE5281",
                        "data/GSE5281/GSE5281_series_matrix_ad.csv",
                        "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc_sce_to_bulk, cv.fit,
                        "results/sc_bulk/", "GSE33000",
                        "data/GSE33000/GSE33000_series_matrix_ad.csv",
                        "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc_sce_to_bulk, cv.fit,
                        "results/sc_bulk/", "GSE104704", 
                        "data/GSE104704/GSE104704_RNAseq_ad.csv",
                        "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc_sce_to_bulk, cv.fit,
                        "results/sc_bulk/", "GSE159699", 
                        "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                        "data/GSE159699/GSE159699_RNAseq_count_nc.csv")

# net --------------------------------------------------------------------

RR_inter <- read.table("../iscience/Download_data_RR.txt", header = T, sep = "\t", fill = T, quote = "")
RR_inter <- subset(RR_inter,subset =
                     RR_inter$Species1 == "Homo sapiens"
                   & RR_inter$Species2 == "Homo sapiens"
                   & RR_inter$Interactor1.Symbol %in% ad_gene
                   & RR_inter$Interactor2.Symbol %in% ad_gene)
nrow(RR_inter)
# no links in RR_inter

RP_inter <- read.table("../iscience/Download_data_RP.txt", header = T, sep = "\t", fill = T, quote = "")
RP_inter <- subset(RP_inter,subset =
                     RP_inter$Species1 == "Homo sapiens"
                   & RP_inter$Species2 == "Homo sapiens"
                   & RP_inter$Interactor1.Symbol %in% ad_gene
                   & RP_inter$Interactor2.Symbol %in% ad_gene)
nrow(RP_inter)

RP_inter <- RP_inter[, c(2, 5)]
RP_inter$Interactor2.Symbol <- paste0(RP_inter$Interactor2.Symbol, " ")

# P_P network data downloaded from string
P_P <- read.table("results/string_interactions_short.tsv")
P_P <- P_P[, 1:2]
P_P$V1 <- paste0(P_P$V1, " ")
P_P$V2 <- paste0(P_P$V2, " ")
colnames(P_P) <- colnames(RP_inter)

links <- rbind(RP_inter, P_P)

nodes <- c(c(links$Interactor1.Symbol), c(links$Interactor2.Symbol)) %>% 
  unique() %>% 
  data.frame()
nodes$type <- "mRNA"
nodes$type[grep(" ", nodes$.)] <- "Protein"

network <- graph_from_data_frame(d = links, vertices = nodes, directed = F)
deg <- degree(network, mode = "all")
coul  <- brewer.pal(3, "Set1")
my_color <- coul[as.numeric(as.factor(V(network)$type))]
{
  png("results/sc_bulk/net.png", width = 1500, height = 1500, res = 300)
  plot(network, 
       layout = layout_as_tree(network, circular = TRUE),
       edge.lty = 2,
       vertex.size = 5, 
       vertex.label.cex = 0.5, 
       vertex.label.famliy = "Times", 
       vertex.label.dist = 1,
       vertex.color = my_color)
  legend("bottomright", 
         # inset=.05, 
         # title= NA,
         c("RNA","Protein"),
         # lty=c(1, 2), 
         pch=19,
         bty = "n",
         cex = 0.6,
         col=c("#E41A1C", "#377EB8"))
  dev.off()
}
netdata <- data.frame(ad_gene)

netdata$links <- 0
for(i in 1:nrow(netdata)) {
  n <- grep(paste0("^", netdata$ad_gene[i]), c(c(links$Interactor1.Symbol), c(links$Interactor2.Symbol))) %>% 
    length()
  netdata$links[i] <- n
}
netdata <- which(netdata$links != 0) %>% netdata[., ]

netdata <- netdata[order(netdata$links, decreasing = F), ]
write.csv(netdata, "results/sc_bulk/netdata_sce.csv", row.names = F)

netdata$ad_gene <- factor(netdata$ad_gene, levels = netdata$ad_gene)
{
  ggplot(netdata, aes(y = ad_gene, x = links, fill = -links)) +
    geom_col() +
    ylab(label = "Gene") +
    xlab(label = "The number of edges") +
    guides(fill = "none") +
    scale_fill_distiller(palette = "Reds") +
    theme_classic()
  ggsave("results/sc_bulk/gene_net_freq.png", width = 3, height = 3, device = "png")
}

# 染色体 ---------------------------------------------------------------------

{
  pdf("results/sc_bulk/chome_gene.pdf", width = 2400, height = 2400)
  
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$grid <- 1
  RCircos.Reset.Plot.Parameters(rcircos.params)
  RCircos.List.Plot.Parameters()
  data("UCSC.HG38.Human.CytoBandIdeogram")
  cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
  RCircos.Set.Core.Components(cyto.info,
                              chr.exclude = NULL, 
                              tracks.inside = 5,
                              tracks.outside = 1)
  
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  
  data("RCircos.Gene.Label.Data")
  name.col <- 4
  side <- "in"
  track.num <- 1
  
  
  
  ELAVL3 <- c("chr19", 11451326, 11481046, "ELAVL3")
  GFAP <- c("chr17", 44903159, 44915500, "GFAP")
  NEUROD6 <- c("chr7", 31337465, 31340726, "NEUROD6")
  RTN1 <- c("chr14", 59595976, 59870776, "RTN1")
  SULT4A1 <- c("chr22", 43824509, 43862513, "SULT4A1")
  FAIM2 <- c("chr12", 49866896, 49903900, "FAIM2")
  RPH3A <- c("chr12", 112575236, 112898881, "RPH3A")
  SLC24A2 <- c("chr9", 19507455, 20307892, "SLC24A2")
  HMP19 <- c("chr5", 174045706, 174109179, "HMP19")
  SNAP25 <- c("chr20", 10218830, 10307418, "SNAP25")
  LRRTM2 <- c("chr5", 138868921, 138875335, "LRRTM2")
  
  RCircos.Gene.Label.Data_1 <- data.frame(LRRTM2,
                                          ELAVL3,
                                          GFAP,
                                          NEUROD6,
                                          RTN1,
                                          SULT4A1,
                                          FAIM2,
                                          RPH3A,
                                          SLC24A2,
                                          HMP19 ,
                                          SNAP25) %>%
    t() %>% 
    data.frame()
  
  colnames(RCircos.Gene.Label.Data_1) <- colnames(RCircos.Gene.Label.Data)
  rownames(RCircos.Gene.Label.Data_1) <- 1:nrow(RCircos.Gene.Label.Data_1)
  RCircos.Gene.Label.Data_1$chromStart <- as.numeric(RCircos.Gene.Label.Data_1$chromStart)
  RCircos.Gene.Label.Data_1$chromEnd <- as.numeric(RCircos.Gene.Label.Data_1$chromEnd)
  # RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data_1, 1, "out")
  RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data_1, 4, 1, "out")
  
  
  data("RCircos.Link.Data")
  head(RCircos.Link.Data)
  
  RCircos.Link.Data_1 <- data.frame(c(ELAVL3, FAIM2),
                                    c(GFAP, SULT4A1),
                                    c(ELAVL3, RPH3A),
                                    c(RTN1, SNAP25),
                                    c(LRRTM2, HMP19),
                                    c(NEUROD6, SLC24A2)) %>%
    t() %>% 
    data.frame()
  RCircos.Link.Data_1 <- RCircos.Link.Data_1[, -c(4, 8)]
  colnames(RCircos.Link.Data_1) <- colnames(RCircos.Link.Data)
  rownames(RCircos.Link.Data_1) <- 1:nrow(RCircos.Link.Data_1)
  RCircos.Link.Data_1$chromStart <- as.numeric(RCircos.Link.Data_1$chromStart)
  RCircos.Link.Data_1$chromEnd <- as.numeric(RCircos.Link.Data_1$chromEnd)
  RCircos.Link.Data_1$chromStart.1 <- as.numeric(RCircos.Link.Data_1$chromStart.1)
  RCircos.Link.Data_1$chromEnd.1 <- as.numeric(RCircos.Link.Data_1$chromEnd.1)
  
  RCircos.Link.Plot(RCircos.Link.Data_1, 1, T)
  dev.off()
  
}

# GO and KEGG----------------------------------------------------------------------

geneid <- bitr(ad_gene, 
               fromType = "SYMBOL",
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)

all_enrich <- enrichGO(geneid$ENTREZID,
                       keyType = "ENTREZID",
                       OrgDb = org.Hs.eg.db, 
                       ont = "ALL",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
all_enrich <- setReadable(all_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
png("results/sc_bulk/go-enrich-net.png", width = 2700, height = 2100, res = 300)

enrichplot::cnetplot(all_enrich,
                     circular = FALSE,
                     showCategory = 10,
                     layout = "circle",
                     # node_label="gene",
                     colorEdge = TRUE)
dev.off()

png("results/sc_bulk/all.png", width = 2700, height = 3600, res = 300)

barplot(all_enrich, drop = TRUE, showCategory = 10, split = "ONTOLOGY", ) + 
  facet_grid(ONTOLOGY~., scale = 'free')

dev.off()
library(clusterProfiler)
cnetplot(all_enrich, categorySize = "pvalue", node_label = "all")
{
  bp_enrich <- enrichGO(geneid$ENTREZID,
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1)
  png("results/sc_bulk/bp.png", width = 300, height = 300)
  dotplot(bp_enrich, title = "BP", )
  dev.off()
  cc_enrich <- enrichGO(geneid$ENTREZID,
                        OrgDb = org.Hs.eg.db, 
                        ont = "CC",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1)
  png("results/sc_bulk/cc.png", width = 1800, height = 1800, res = 300)
  dotplot(cc_enrich, title = "CC")
  dev.off()
  mf_enrich <- enrichGO(geneid$ENTREZID,
                        OrgDb = org.Hs.eg.db, 
                        ont = "MF",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1)
  png("results/sc_bulk/mf.png", width = 1800, height = 1800, res = 300)
  dotplot(mf_enrich, title = "MF")
  dev.off()
}
# both of enrich have no results
# don't need to write functions

# expression level --------------------------------------------------------

lasso_pair_slc_sce_to_bulk[5, ] <- lasso_pair_slc_sce_to_bulk[5, 2:1]
# percent bar plot
xxs_exp_col_plot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE48350", 
                 "data/GSE48350/GSE48350_series_matrix_ad.csv",
                 "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_exp_col_plot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE5281",
                 "data/GSE5281/GSE5281_series_matrix_ad.csv",
                 "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_exp_col_plot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE33000",
                 "data/GSE33000/GSE33000_series_matrix_ad.csv",
                 "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_exp_col_plot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE104704", 
                 "data/GSE104704/GSE104704_RNAseq_ad.csv",
                 "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_exp_col_plot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE159699", 
                 "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                 "data/GSE159699/GSE159699_RNAseq_count_nc.csv")

# waterfall --------------------------------------------------------------

xxs_waterfall_6(lasso_pair_slc_sce_to_bulk, "GSE48350", "results/sc_bulk/",
              "data/GSE48350/GSE48350_series_matrix_ad.csv",
              "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_waterfall_6(lasso_pair_slc_sce_to_bulk, "GSE5281", "results/sc_bulk/",
              "data/GSE5281/GSE5281_series_matrix_ad.csv",
              "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_waterfall_6(lasso_pair_slc_sce_to_bulk, "GSE33000", "results/sc_bulk/",
              "data/GSE33000/GSE33000_series_matrix_ad.csv",
              "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_waterfall_6(lasso_pair_slc_sce_to_bulk, "GSE104704", "results/sc_bulk/",
              "data/GSE104704/GSE104704_RNAseq_ad.csv",
              "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_waterfall_6(lasso_pair_slc_sce_to_bulk, "GSE159699", "results/sc_bulk/",
              "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
              "data/GSE159699/GSE159699_RNAseq_count_nc.csv")

# boxplot ---------------------------------------------------------------

xxs_boxplot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE48350", 
            "data/GSE48350/GSE48350_series_matrix_ad.csv",
            "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_boxplot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE5281",
            "data/GSE5281/GSE5281_series_matrix_ad.csv",
            "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_boxplot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE33000",
            "data/GSE33000/GSE33000_series_matrix_ad.csv",
            "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_boxplot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE104704", 
            "data/GSE104704/GSE104704_RNAseq_ad.csv",
            "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_boxplot(lasso_pair_slc_sce_to_bulk, "results/sc_bulk/", "GSE159699", 
            "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
            "data/GSE159699/GSE159699_RNAseq_count_nc.csv")

# umap -----------------------------------------------------------------

umap_location <- combined_sce@reductions[["umap"]]@cell.embeddings %>% # 这个是拿到所有细胞
  data.frame()
cell_annotation <- combined_sce@meta.data[["celltype"]] %>% 
  data.frame()
umap_location <- cbind(umap_location, cell_annotation)

ad_location <- grep("_1", rownames(umap_location)) %>% 
  umap_location[., ]
nc_location <- grep("_2", rownames(umap_location)) %>% 
  umap_location[., ]

ad_adj_matrix <- xxs_get_adj_matrix_1(lasso_pair_slc_sce_to_bulk, sce_ad_seurat@assays$RNA@counts %>%
                                        as.matrix()) %>% 
  data.frame()
nc_adj_matrix <- xxs_get_adj_matrix_1(lasso_pair_slc_sce_to_bulk, sce_nc_seurat@assays$RNA@counts %>%
                                        as.matrix()) %>% 
  data.frame()

ad_adj_matrix <- ifelse(ad_adj_matrix == 1, "Positive", "Negative")
nc_adj_matrix <- ifelse(nc_adj_matrix == 1, "Positive", "Negative")
ad_plot_data <- cbind(ad_location, ad_adj_matrix)
for(i in 4:9) {
  ad_plot_data[, i] <- factor(ad_plot_data[, i], levels = c("Positive", "Negative"))
}

nc_plot_data <- cbind(nc_location, nc_adj_matrix)
for(i in 4:9) {
  nc_plot_data[, i] <- factor(nc_plot_data[, i], levels = c("Positive", "Negative"))
}

# Umap 
{
  AD_title <- round(sum(ad_plot_data$ELAVL3_vs_FAIM2 == "Positive") / nrow(ad_plot_data) * 100, 2) %>% 
    paste0(., "% in AD")
  NC_title <- round(sum(nc_plot_data$ELAVL3_vs_FAIM2 == "Positive") / nrow(ad_plot_data) * 100, 2) %>% 
    paste0(., "% in NC")
  p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = ELAVL3_vs_FAIM2), size = 0.1) +
    scale_colour_manual(values = c("red3", "gray")) +
    theme_classic() +
    guides(color = "none") +
    ggtitle(AD_title)
  p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = ELAVL3_vs_FAIM2), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    ggtitle(NC_title)
  p <- p1 + p2
  
  ggsave(filename = "results/sc_bulk/ELAVL3_vs_FAIM2.png", plot = p, device = "png", dpi = 300, width = 6, height = 2.5)
}
{
  AD_title <- round(sum(ad_plot_data$ELAVL3_vs_RPH3A == "Positive") / nrow(ad_plot_data) * 100, 2) %>% 
    paste0(., "% in AD")
  NC_title <- round(sum(nc_plot_data$ELAVL3_vs_RPH3A == "Positive") / nrow(ad_plot_data) * 100, 2) %>% 
    paste0(., "% in NC")
  p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = ELAVL3_vs_RPH3A), size = 0.1) +
    theme_classic() +
    guides(color = "none") +
    scale_colour_manual(values = c("red3", "gray")) +
    guides(color = "none") +
    ggtitle(AD_title)
  p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = ELAVL3_vs_RPH3A), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    ggtitle(NC_title)
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/ELAVL3_vs_RPH3A.png", plot = p, device = "png", dpi = 300, width = 6, height = 2.5)
}
{
  p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = GFAP_vs_SULT4A1), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    guides(color = "none") +
    ggtitle("AD")
  
  p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = GFAP_vs_SULT4A1), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    ggtitle("NC")
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/GFAP_vs_SULT4A1.png", plot = p, device = "png", dpi = 300, width = 6, height = 2.5)
}
{
  p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = LRRTM2_vs_HMP19), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    guides(color = "none") +
    ggtitle("AD")
  
  p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = LRRTM2_vs_HMP19), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    ggtitle("NC")
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/LRRTM2_vs_HMP19.png", plot = p, device = "png", dpi = 300, width = 6, height = 2.5)
}
{
  p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = SLC24A2_vs_NEUROD6), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    guides(color = "none") +
    ggtitle("AD")
  
  p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = SLC24A2_vs_NEUROD6), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    ggtitle("NC")
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/NEUROD6_vs_SLC24A2.png", plot = p, device = "png", dpi = 300, width = 6, height = 2.5)
}
{
  AD_title <- round(sum(ad_plot_data$RTN1_vs_SNAP25 == "Positive") / nrow(ad_plot_data) * 100, 2) %>% 
    paste0(., "% in AD")
  NC_title <- round(sum(nc_plot_data$RTN1_vs_SNAP25 == "Positive") / nrow(ad_plot_data) * 100, 2) %>% 
    paste0(., "% in NC")
  p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = RTN1_vs_SNAP25), size = 0.1) +
    theme_classic() +
    guides(color = "none") +
    scale_colour_manual(values = c("red3", "gray")) +
    ggtitle(AD_title)
  p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = RTN1_vs_SNAP25), size = 0.1) +
    theme_classic() +
    scale_colour_manual(values = c("red3", "gray")) +
    ggtitle(NC_title)
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/RTN1_vs_SNAP25.png", plot = p, device = "png", dpi = 300, width = 6, height = 2.5)
  
}

# 饼图
{
  
  
  {
    
    ggplot(ad_plot_data, aes(x = "", y = ELAVL3_vs_FAIM2, fill = ELAVL3_vs_FAIM2)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/ELAVL3_vs_FAIM2_pie_ad.png", device = "png", dpi = 300, width = 3, height = 3)
    
    ggplot(nc_plot_data, aes(x = "", y = ELAVL3_vs_FAIM2, fill = ELAVL3_vs_FAIM2)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/ELAVL3_vs_FAIM2_pie_nc.png", device = "png", dpi = 300, width = 3, height = 3)
    
  }
  {
    ggplot(ad_plot_data, aes(x = "", y = ELAVL3_vs_RPH3A, fill = ELAVL3_vs_RPH3A)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/ELAVL3_vs_RPH3A_pie_ad.png", device = "png", dpi = 300, width = 3, height = 3)
    
    ggplot(nc_plot_data, aes(x = "", y = ELAVL3_vs_RPH3A, fill = ELAVL3_vs_RPH3A)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/ELAVL3_vs_RPH3A_pie_nc.png", device = "png", dpi = 300, width = 3, height = 3)
    
  }
  {
    ggplot(ad_plot_data, aes(x = "", y = GFAP_vs_SULT4A1, fill = GFAP_vs_SULT4A1)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/GFAP_vs_SULT4A1_pie_ad.png", device = "png", dpi = 300, width = 3, height = 3)
    
    ggplot(nc_plot_data, aes(x = "", y = GFAP_vs_SULT4A1, fill = GFAP_vs_SULT4A1)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/GFAP_vs_SULT4A1_pie_nc.png", device = "png", dpi = 300, width = 3, height = 3)
    
  }
  {
    ggplot(ad_plot_data, aes(x = "", y = LRRTM2_vs_HMP19, fill = LRRTM2_vs_HMP19)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/LRRTM2_vs_HMP19_pie_ad.png", device = "png", dpi = 300, width = 3, height = 3)
    
    ggplot(nc_plot_data, aes(x = "", y = LRRTM2_vs_HMP19, fill = LRRTM2_vs_HMP19)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/LRRTM2_vs_HMP19_pie_nc.png", device = "png", dpi = 300, width = 3, height = 3)
    
  }
  {
    ggplot(ad_plot_data, aes(x = "", y = SLC24A2_vs_NEUROD6, fill = SLC24A2_vs_NEUROD6)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/SLC24A2_vs_NEUROD6_pie_ad.png", device = "png", dpi = 300, width = 3, height = 3)
    
    ggplot(nc_plot_data, aes(x = "", y = SLC24A2_vs_NEUROD6, fill = SLC24A2_vs_NEUROD6)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/SLC24A2_vs_NEUROD6_pie_nc.png", device = "png", dpi = 300, width = 3, height = 3)
    
  }
  {
    ggplot(ad_plot_data, aes(x = "", y = RTN1_vs_SNAP25, fill = RTN1_vs_SNAP25)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/RTN1_vs_SNAP25_pie_ad.png", device = "png", dpi = 300, width = 3, height = 3)
    
    ggplot(nc_plot_data, aes(x = "", y = RTN1_vs_SNAP25, fill = RTN1_vs_SNAP25)) +
      geom_bar(stat = 'identity') +
      coord_polar(theta = "y") +
      guides(fill = "none") +
      scale_fill_manual(values = c("red3", "gray")) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
      theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
    ggsave(filename = "results/sc_bulk/RTN1_vs_SNAP25_pie_nc.png", device = "png", dpi = 300, width = 3, height = 3)
    
  }
}

# adj barplot
# ELAVL3_vs_FAIM2
{
  adj_plot_data_ad <- ad_plot_data[, 3:9]
  adj_plot_data_ad <- table(adj_plot_data_ad$., adj_plot_data_ad$ELAVL3_vs_FAIM2) %>% 
    data.frame()
  adj_plot_data_ad$Var2 <- adj_plot_data_ad$Freq[1:4] + adj_plot_data_ad$Freq[5:8]
  adj_plot_data_ad$`Percent (%)` <- adj_plot_data_ad$Freq / adj_plot_data_ad$Var2 * 100
  adj_plot_data_ad$Group <- "AD"
  
  adj_plot_data_nc <- nc_plot_data[, 3:9]
  adj_plot_data_nc <- table(adj_plot_data_nc$., adj_plot_data_nc$ELAVL3_vs_FAIM2) %>% 
    data.frame()
  adj_plot_data_nc$Var2 <- adj_plot_data_nc$Freq[1:4] + adj_plot_data_nc$Freq[5:8]
  adj_plot_data_nc$`Percent (%)` <- adj_plot_data_nc$Freq / adj_plot_data_nc$Var2 * 100
  adj_plot_data_nc$Group <- "NC"
  
  Fisher <- rbind(adj_plot_data_ad[, c(1, 3)], adj_plot_data_nc[, c(1, 3)])
  Fisher_adj <- matrix(Fisher$Freq, nrow = 4, byrow = F) %>% 
    data.frame()
  rownames(Fisher_adj) <- Fisher$Var1[1:4]
  colnames(Fisher_adj) <- c("P_AD", "N_AD", "P_NC", "N_CA")
  Fisher_adj$p_value <- 0
  for (i in 1:4) {
    Fisher_adj$p_value[i] <- 
      Fisher_adj[i, 1:4] %>% 
      as.numeric() %>% 
      matrix(., ncol = 2, nrow = 2) %>% 
      fisher.test() %$%
      p.value
  }
  write.csv(Fisher_adj, "results/sc_bulk/ELAVL3_vs_FAIM2_fisher.csv")
  
  adj_plot_data <- rbind(adj_plot_data_ad[1:4, c(1, 4, 5)], adj_plot_data_nc[1:4,  c(1, 4, 5)])
  colnames(adj_plot_data)[1] <- "Cell Type"
  ggplot(adj_plot_data) +
    geom_bar(aes(y = `Cell Type`, x = `Percent (%)`, fill = Group), stat="identity", position = 'dodge') +
    theme_classic() +
    ggtitle("ELAVL3_vs_FAIM2")
  ggsave("results/sc_bulk/ELAVL3_vs_FAIM2_bar.png", width = 5, height = 3, dpi = 300, device = "png")
  
}

# ELAVL3_vs_RPH3A
{
  adj_plot_data_ad <- ad_plot_data[, 3:9]
  adj_plot_data_ad <- table(adj_plot_data_ad$., adj_plot_data_ad$ELAVL3_vs_RPH3A) %>% 
    data.frame()
  adj_plot_data_ad$Var2 <- adj_plot_data_ad$Freq[1:4] + adj_plot_data_ad$Freq[5:8]
  adj_plot_data_ad$`Percent (%)` <- adj_plot_data_ad$Freq / adj_plot_data_ad$Var2 * 100
  adj_plot_data_ad$Group <- "AD"
  
  adj_plot_data_nc <- nc_plot_data[, 3:9]
  adj_plot_data_nc <- table(adj_plot_data_nc$., adj_plot_data_nc$ELAVL3_vs_RPH3A) %>% 
    data.frame()
  adj_plot_data_nc$Var2 <- adj_plot_data_nc$Freq[1:4] + adj_plot_data_nc$Freq[5:8]
  adj_plot_data_nc$`Percent (%)` <- adj_plot_data_nc$Freq / adj_plot_data_nc$Var2 * 100
  adj_plot_data_nc$Group <- "NC"
  
  Fisher <- rbind(adj_plot_data_ad[, c(1, 3)], adj_plot_data_nc[, c(1, 3)])
  Fisher_adj <- matrix(Fisher$Freq, nrow = 4, byrow = F) %>% 
    data.frame()
  rownames(Fisher_adj) <- Fisher$Var1[1:4]
  colnames(Fisher_adj) <- c("P_AD", "N_AD", "P_NC", "N_CA")
  Fisher_adj$p_value <- 0
  for (i in 1:4) {
    Fisher_adj$p_value[i] <- 
      Fisher_adj[i, 1:4] %>% 
      as.numeric() %>% 
      matrix(., ncol = 2, nrow = 2) %>% 
      fisher.test() %$%
      p.value
  }
  write.csv(Fisher_adj, "results/sc_bulk/ELAVL3_vs_RPH3A_fisher.csv")
  
  adj_plot_data <- rbind(adj_plot_data_ad[1:4, c(1, 4, 5)], adj_plot_data_nc[1:4,  c(1, 4, 5)])
  colnames(adj_plot_data)[1] <- "Cell Type"
  ggplot(adj_plot_data) +
    geom_bar(aes(y = `Cell Type`, x = `Percent (%)`, fill = Group), stat="identity", position = 'dodge') +
    theme_classic() +
    ggtitle("ELAVL3_vs_RPH3A")
  ggsave("results/sc_bulk/ELAVL3_vs_RPH3A_bar.png", width = 5, height = 3, dpi = 300, device = "png")
}

# GFAP_vs_SULT4A1
{
  adj_plot_data_ad <- ad_plot_data[, 3:9]
  adj_plot_data_ad <- table(adj_plot_data_ad$., adj_plot_data_ad$GFAP_vs_SULT4A1) %>% 
    data.frame()
  adj_plot_data_ad$Var2 <- adj_plot_data_ad$Freq[1:4] + adj_plot_data_ad$Freq[5:8]
  adj_plot_data_ad$`Percent (%)` <- adj_plot_data_ad$Freq / adj_plot_data_ad$Var2 * 100
  adj_plot_data_ad$Group <- "AD"
  
  adj_plot_data_nc <- nc_plot_data[, 3:9]
  adj_plot_data_nc <- table(adj_plot_data_nc$., adj_plot_data_nc$GFAP_vs_SULT4A1) %>% 
    data.frame()
  adj_plot_data_nc$Var2 <- adj_plot_data_nc$Freq[1:4] + adj_plot_data_nc$Freq[5:8]
  adj_plot_data_nc$`Percent (%)` <- adj_plot_data_nc$Freq / adj_plot_data_nc$Var2 * 100
  adj_plot_data_nc$Group <- "NC"
  
  Fisher <- rbind(adj_plot_data_ad[, c(1, 3)], adj_plot_data_nc[, c(1, 3)])
  Fisher_adj <- matrix(Fisher$Freq, nrow = 4, byrow = F) %>% 
    data.frame()
  rownames(Fisher_adj) <- Fisher$Var1[1:4]
  colnames(Fisher_adj) <- c("P_AD", "N_AD", "P_NC", "N_CA")
  Fisher_adj$p_value <- 0
  for (i in 1:4) {
    Fisher_adj$p_value[i] <- 
      Fisher_adj[i, 1:4] %>% 
      as.numeric() %>% 
      matrix(., ncol = 2, nrow = 2) %>% 
      fisher.test() %$%
      p.value
  }
  write.csv(Fisher_adj, "results/sc_bulk/GFAP_vs_SULT4A1_fisher.csv")
  
  adj_plot_data <- rbind(adj_plot_data_ad[1:4, c(1, 4, 5)], adj_plot_data_nc[1:4,  c(1, 4, 5)])
  colnames(adj_plot_data)[1] <- "Cell Type"
  ggplot(adj_plot_data) +
    geom_bar(aes(y = `Cell Type`, x = `Percent (%)`, fill = Group), stat="identity", position = 'dodge') +
    theme_classic() +
    ggtitle("GFAP_vs_SULT4A1")
  ggsave("results/sc_bulk/GFAP_vs_SULT4A1_bar.png", width = 5, height = 3, dpi = 300, device = "png")
}

# LRRTM2_vs_HMP19
{
  adj_plot_data_ad <- ad_plot_data[, 3:9]
  adj_plot_data_ad <- table(adj_plot_data_ad$., adj_plot_data_ad$LRRTM2_vs_HMP19) %>% 
    data.frame()
  adj_plot_data_ad$Var2 <- adj_plot_data_ad$Freq[1:4] + adj_plot_data_ad$Freq[5:8]
  adj_plot_data_ad$`Percent (%)` <- adj_plot_data_ad$Freq / adj_plot_data_ad$Var2 * 100
  adj_plot_data_ad$Group <- "AD"
  
  adj_plot_data_nc <- nc_plot_data[, 3:9]
  adj_plot_data_nc <- table(adj_plot_data_nc$., adj_plot_data_nc$LRRTM2_vs_HMP19) %>% 
    data.frame()
  adj_plot_data_nc$Var2 <- adj_plot_data_nc$Freq[1:4] + adj_plot_data_nc$Freq[5:8]
  adj_plot_data_nc$`Percent (%)` <- adj_plot_data_nc$Freq / adj_plot_data_nc$Var2 * 100
  adj_plot_data_nc$Group <- "NC"
  
  Fisher <- rbind(adj_plot_data_ad[, c(1, 3)], adj_plot_data_nc[, c(1, 3)])
  Fisher_adj <- matrix(Fisher$Freq, nrow = 4, byrow = F) %>% 
    data.frame()
  rownames(Fisher_adj) <- Fisher$Var1[1:4]
  colnames(Fisher_adj) <- c("P_AD", "N_AD", "P_NC", "N_CA")
  Fisher_adj$p_value <- 0
  for (i in 1:4) {
    Fisher_adj$p_value[i] <- 
      Fisher_adj[i, 1:4] %>% 
      as.numeric() %>% 
      matrix(., ncol = 2, nrow = 2) %>% 
      fisher.test() %$%
      p.value
  }
  write.csv(Fisher_adj, "results/sc_bulk/LRRTM2_vs_HMP19_fisher.csv")
  
  adj_plot_data <- rbind(adj_plot_data_ad[1:4, c(1, 4, 5)], adj_plot_data_nc[1:4,  c(1, 4, 5)])
  colnames(adj_plot_data)[1] <- "Cell Type"
  ggplot(adj_plot_data) +
    geom_bar(aes(y = `Cell Type`, x = `Percent (%)`, fill = Group), stat="identity", position = 'dodge') +
    theme_classic() +
    ggtitle("LRRTM2_vs_HMP19")
  ggsave("results/sc_bulk/LRRTM2_vs_HMP19_bar.png", width = 5, height = 3, dpi = 300, device = "png")
}

# SLC24A2_vs_NEUROD6
{
  adj_plot_data_ad <- ad_plot_data[, 3:9]
  adj_plot_data_ad <- table(adj_plot_data_ad$., adj_plot_data_ad$SLC24A2_vs_NEUROD6) %>% 
    data.frame()
  adj_plot_data_ad$Var2 <- adj_plot_data_ad$Freq[1:4] + adj_plot_data_ad$Freq[5:8]
  adj_plot_data_ad$`Percent (%)` <- adj_plot_data_ad$Freq / adj_plot_data_ad$Var2 * 100
  adj_plot_data_ad$Group <- "AD"
  
  adj_plot_data_nc <- nc_plot_data[, 3:9]
  adj_plot_data_nc <- table(adj_plot_data_nc$., adj_plot_data_nc$SLC24A2_vs_NEUROD6) %>% 
    data.frame()
  adj_plot_data_nc$Var2 <- adj_plot_data_nc$Freq[1:4] + adj_plot_data_nc$Freq[5:8]
  adj_plot_data_nc$`Percent (%)` <- adj_plot_data_nc$Freq / adj_plot_data_nc$Var2 * 100
  adj_plot_data_nc$Group <- "NC"
  
  Fisher <- rbind(adj_plot_data_ad[, c(1, 3)], adj_plot_data_nc[, c(1, 3)])
  Fisher_adj <- matrix(Fisher$Freq, nrow = 4, byrow = F) %>% 
    data.frame()
  rownames(Fisher_adj) <- Fisher$Var1[1:4]
  colnames(Fisher_adj) <- c("P_AD", "N_AD", "P_NC", "N_CA")
  Fisher_adj$p_value <- 0
  for (i in 1:4) {
    Fisher_adj$p_value[i] <- 
      Fisher_adj[i, 1:4] %>% 
      as.numeric() %>% 
      matrix(., ncol = 2, nrow = 2) %>% 
      fisher.test() %$%
      p.value
  }
  write.csv(Fisher_adj, "results/sc_bulk/SLC24A2_vs_NEUROD6_fisher.csv")
  
  adj_plot_data <- rbind(adj_plot_data_ad[1:4, c(1, 4, 5)], adj_plot_data_nc[1:4,  c(1, 4, 5)])
  colnames(adj_plot_data)[1] <- "Cell Type"
  ggplot(adj_plot_data) +
    geom_bar(aes(y = `Cell Type`, x = `Percent (%)`, fill = Group), stat="identity", position = 'dodge') +
    theme_classic() +
    ggtitle("SLC24A2_vs_NEUROD6")
  ggsave("results/sc_bulk/SLC24A2_vs_NEUROD6_bar.png", width = 5, height = 3, dpi = 300, device = "png")
}

# RTN1_vs_SNAP25
{
  adj_plot_data_ad <- ad_plot_data[, 3:9]
  adj_plot_data_ad <- table(adj_plot_data_ad$., adj_plot_data_ad$RTN1_vs_SNAP25) %>% 
    data.frame()
  adj_plot_data_ad$Var2 <- adj_plot_data_ad$Freq[1:4] + adj_plot_data_ad$Freq[5:8]
  adj_plot_data_ad$`Percent (%)` <- adj_plot_data_ad$Freq / adj_plot_data_ad$Var2 * 100
  adj_plot_data_ad$Group <- "AD"
  
  adj_plot_data_nc <- nc_plot_data[, 3:9]
  adj_plot_data_nc <- table(adj_plot_data_nc$., adj_plot_data_nc$RTN1_vs_SNAP25) %>% 
    data.frame()
  adj_plot_data_nc$Var2 <- adj_plot_data_nc$Freq[1:4] + adj_plot_data_nc$Freq[5:8]
  adj_plot_data_nc$`Percent (%)` <- adj_plot_data_nc$Freq / adj_plot_data_nc$Var2 * 100
  adj_plot_data_nc$Group <- "NC"
  
  Fisher <- rbind(adj_plot_data_ad[, c(1, 3)], adj_plot_data_nc[, c(1, 3)])
  Fisher_adj <- matrix(Fisher$Freq, nrow = 4, byrow = F) %>% 
    data.frame()
  rownames(Fisher_adj) <- Fisher$Var1[1:4]
  colnames(Fisher_adj) <- c("P_AD", "N_AD", "P_NC", "N_CA")
  Fisher_adj$p_value <- 0
  for (i in 1:4) {
    Fisher_adj$p_value[i] <- 
      Fisher_adj[i, 1:4] %>% 
      as.numeric() %>% 
      matrix(., ncol = 2, nrow = 2) %>% 
      fisher.test() %$%
      p.value
  }
  write.csv(Fisher_adj, "results/sc_bulk/RTN1_vs_SNAP25_fisher.csv")
  
  adj_plot_data <- rbind(adj_plot_data_ad[1:4, c(1, 4, 5)], adj_plot_data_nc[1:4,  c(1, 4, 5)])
  colnames(adj_plot_data)[1] <- "Cell Type"
  ggplot(adj_plot_data) +
    geom_bar(aes(y = `Cell Type`, x = `Percent (%)`, fill = Group), stat="identity", position = 'dodge') +
    theme_classic() +
    ggtitle("RTN1_vs_SNAP25")
  ggsave("results/sc_bulk/RTN1_vs_SNAP25_bar.png", width = 5, height = 3, dpi = 300, device = "png")
  
  ggplot(adj_plot_data) +
    geom_bar(aes(y = `Cell Type`, x = `Percent (%)`, fill = Group), stat="identity", position = 'dodge') +
    theme_classic() +
    ggtitle("RTN1_vs_SNAP25")
  ggsave("results/sc_bulk/RTN1_vs_SNAP25_bar-label.png", width = 5, height = 3, dpi = 300, device = "png")
  
}

