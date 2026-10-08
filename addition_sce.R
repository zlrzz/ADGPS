# 加载包 ---------------------------------------------------------------------

source("final-code/analysis_data_functions.R", encoding = "UTF-8")

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


astrocytes <- c("GFAP", "AQP4", "SLC1A2")
excitatory <- c("CAMK2A", "LDB2", "SLC17A7")
inhibitory <- c("GAD1", "GAD2")
oligodendrocytes <- c("MBP", "MOG", "PLP1")

gene_anno_list <- list(excitatory = excitatory,
                       inhibitory = inhibitory,
                       astrocytes = astrocytes,
                       oligodendrocytes = oligodendrocytes)

Idents(combined_sce) = combined_sce$integrated_snn_res.1
png("results/sc/DotPlot_gene_anno_list.png", width = 4500, height = 1500, res = 300)
DotPlot(combined_sce, features = gene_anno_list, assay = "RNA")
dev.off()
  

Astro <- 21
Oligo <- 20
Inhib <- c(0, 6, 10, 14, 16, 19)
Excit <- setdiff(c(1:19), Inhib)

combined_sce$seurat_clusters <- factor(combined_sce$seurat_clusters, levels = c(Astro, Oligo, Excit, Inhib))
combined_sce$celltype = ifelse(combined_sce$seurat_clusters %in% Astro, "astrocytes",
                               ifelse(combined_sce$seurat_clusters %in% Oligo, "oligodendrocytes",
                                      ifelse(combined_sce$seurat_clusters %in% Excit, "excitatory",
                                             "inhibitory")))
Idents(combined_sce) = combined_sce$celltype
combined_sce$celltype <- factor(combined_sce$celltype, levels = c("excitatory", "inhibitory", "astrocytes", "oligodendrocytes"))


png("results/sc/VlnPlot.gene_anno.png", width = 2100, height = 3200, res = 300)
VlnPlot(combined_sce, features = unlist(gene_anno_list), group.by = "celltype", assay = "RNA", pt.size = 0, ncol = 3)
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

set.seed(38)
train_index_ad <- sample(1:ncol(sce_ad),  ncol(sce_ad) * 0.7)
train_index_nc <- sample(1:ncol(sce_nc),  ncol(sce_nc) * 0.7)
sce_ad_test <- sce_ad[, -train_index_ad]
sce_nc_test <- sce_nc[, -train_index_nc]
sce_ad_train <- sce_ad[, train_index_ad]
sce_nc_train <- sce_nc[, train_index_nc]

# 内部拆分为训练和验证组 -------------------------------------------------------------

pair_data_sce <- xxs_get_pair_data_from_matrix(sce_ad_train, sce_nc_train, rate_cutoff = 0.6)
lasso_pair_sce <- xxs_get_pair_from_lasso(pair_data_sce, sce_ad_train, sce_nc_train)

lasso_pair_plot_sce <- subset(lasso_pair_sce, subset = lasso_pair_sce$Freq > 30)
lasso_pair_plot_sce <- lasso_pair_plot_sce[order(lasso_pair_plot_sce$Freq, decreasing = F), ]
lasso_pair_plot_sce$. <- factor(lasso_pair_plot_sce$., levels = lasso_pair_plot_sce$.)
colnames(lasso_pair_plot_sce) <- c("Gene Pair", "Freq")
{
  ggplot(lasso_pair_plot_sce) +
    geom_point(aes(x = Freq, y = `Gene Pair`, col = `Gene Pair`)) +
    guides(col = "none") +
    theme_classic() +
    theme(axis.text.y = element_blank())
  
  ggsave("results/sc/gene_pair_freq.png", device = "png", width = 4, height = 5, dpi = 300)
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
Active.Coefficients_sce$`Gene Pair` <- factor(row.names(Active.Coefficients_sce), levels = row.names(Active.Coefficients_sce)[rank])
print(Active.Coefficients_sce)
{
  ggplot(Active.Coefficients_sce) +
    geom_point(aes(x = Active.Coefficients_sce, y = `Gene Pair`, col = Active.Coefficients_sce)) +
    guides(col = "none") +
    theme_classic() +
    scale_color_gradient(low = "blue", high = "red") +
    xlab("Active.Coefficients") +
    theme(axis.text.y = element_blank())
  
  ggsave("results/sc/slc_pair_coeff.png", device = "png", width = 4, height = 5, dpi = 300)
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
pair_data_from_sce$gene_pair <- Active.Coefficients_sce$`Gene Pair`

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

lasso_pair_plot_sce_to_bulk <- subset(lasso_pair_bulk_sce_to_bulk, subset = lasso_pair_bulk_sce_to_bulk$Freq > 30)
lasso_pair_plot_sce_to_bulk <- lasso_pair_plot_sce_to_bulk[order(lasso_pair_plot_sce_to_bulk$Freq, decreasing = F), ]
lasso_pair_plot_sce_to_bulk$. <- factor(lasso_pair_plot_sce_to_bulk$., levels = lasso_pair_plot_sce_to_bulk$.)
colnames(lasso_pair_plot_sce_to_bulk) <- c("Gene Pair", "Freq")

{
  ggplot(lasso_pair_plot_sce_to_bulk) +
    geom_col(aes(x = Freq, y = `Gene Pair`, fill = `Gene Pair`)) +
    guides(fill = "none") +
    theme_classic()
  
  ggsave("results/sc_bulk/gene_pair_freq.png", device = "png", width = 4, height = 5, dpi = 300)
}

lasso_pair_slc_sce_to_bulk <- subset(lasso_pair_bulk_sce_to_bulk, subset = lasso_pair_bulk_sce_to_bulk$Freq > 60, select = .) %>% 
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
rank <- order(Active.Coefficients$Active.Coefficients, decreasing = F)
Active.Coefficients$`Gene Pair` <- factor(row.names(Active.Coefficients), levels = row.names(Active.Coefficients)[rank])
print(Active.Coefficients)
write.table(Active.Coefficients, 
            "results/sc_bulk/Active.Coefficients_train_from_bulk.txt",
            row.names = F)

{ 
  ggplot(Active.Coefficients) +
    geom_col(aes(x = Active.Coefficients, y = `Gene Pair`, fill = -Active.Coefficients)) +
    guides(fill = "none") +
    scale_fill_distiller(palette = "Reds") +
    theme_classic()
  ggsave("results/sc_bulk/slc_pair_coeff.png", device = "png", width = 4, height = 5, dpi = 300)
}

pred_train <- predict(cv.fit, matrix_cbind_train, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
{
  png("results/sc_bulk/bulk pair data train auc.png", width = 900, height = 900, units = "px", res = 300)
  auc1 <- roc(label_train, pred_train, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Training")
  dev.off()
  }# plot(auc1, main = round(auc1$auc, 4))


pred_test <- predict(cv.fit, matrix_cbind_test, type="response", s = Lambda) %>%  #s=seLambda
  as.numeric()
{
  png("results/sc_bulk/bulk pair data test auc.png", width = 900, height = 900, units = "px", res = 300)
  auc2 <- roc(label_test, pred_test, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Testing")
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
RP_inter$Interactor2.Symbol <- paste0(RP_inter$Interactor2.Symbol, "(P)")

# P_P network data downloaded from string
P_P <- read.table("results/string_interactions_short.tsv")
P_P <- P_P[, 1:2]
P_P$V1 <- paste0(P_P$V1, "(P)")
P_P$V2 <- paste0(P_P$V2, "(P)")
colnames(P_P) <- colnames(RP_inter)

links <- rbind(RP_inter, P_P)

nodes <- c(c(links$Interactor1.Symbol), c(links$Interactor2.Symbol)) %>% 
  unique() %>% 
  data.frame()
nodes$type <- "mRNA"
nodes$type[grep("\\(P)", nodes$.)] <- "Protein"

network <- graph_from_data_frame(d = links, vertices = nodes, directed = F)
deg <- degree(network, mode = "all")
coul  <- brewer.pal(3, "Set1")
my_color <- coul[as.numeric(as.factor(V(network)$type))]
{
  png("results/sc_bulk/net.png", width = 1500, height = 1500, res = 300)
  plot(network, 
     vertex.size = 24, 
     vertex.label.cex = 0.7, 
     vertex.label.dist = 3,
     vertex.color = my_color)
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
write.csv(netdata, "results/netdata_sce.csv", row.names = F)

netdata$ad_gene <- factor(netdata$ad_gene, levels = netdata$ad_gene)
{
  ggplot(netdata, aes(y = ad_gene, x = links, fill = -links)) +
    geom_col() +
    ylab(label = "gene") +
    guides(fill = "none") +
    scale_fill_distiller(palette = "Reds") +
    theme_classic()
  ggsave("results/sc_bulk//gene_net_freq.png", width = 3, height = 3, device = "png")
}

# 染色体 ---------------------------------------------------------------------

{
  png("results/sc_bulk/chome_gene.png", width = 1200, height = 1200, res = 300)
  
data("UCSC.HG38.Human.CytoBandIdeogram")
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info,
                            chr.exclude = NULL, 
                            tracks.inside = 2,
                            tracks.outside = 2)

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

RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data_1, 1, "out")
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data_1, 4, 2, "out")


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

bp_enrich <- enrichGO(geneid$ENTREZID,
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)
png("results/sc_bulk/bp.png", width = 1800, height = 1800, res = 300)
dotplot(bp_enrich, title = "BP")
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

xxs_waterfall(lasso_pair_slc_sce_to_bulk, "GSE48350", "results/sc_bulk/",
              "data/GSE48350/GSE48350_series_matrix_ad.csv",
              "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_waterfall(lasso_pair_slc_sce_to_bulk, "GSE5281", "results/sc_bulk/",
              "data/GSE5281/GSE5281_series_matrix_ad.csv",
              "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_waterfall(lasso_pair_slc_sce_to_bulk, "GSE33000", "results/sc_bulk/",
              "data/GSE33000/GSE33000_series_matrix_ad.csv",
              "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_waterfall(lasso_pair_slc_sce_to_bulk, "GSE104704", "results/sc_bulk/",
              "data/GSE104704/GSE104704_RNAseq_ad.csv",
              "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_waterfall(lasso_pair_slc_sce_to_bulk, "GSE159699", "results/sc_bulk/",
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

ad_adj_matrix <- ifelse(ad_adj_matrix == 1, ">", ifelse(ad_adj_matrix == 0, "=", "<"))
nc_adj_matrix <- ifelse(nc_adj_matrix == 1, ">", ifelse(nc_adj_matrix == 0, "=", "<"))
ad_plot_data <- cbind(ad_location, ad_adj_matrix)
for(i in 4:9) {
  ad_plot_data[, i] <- factor(ad_plot_data[, i], levels = c(">", "=", "<"))
}

nc_plot_data <- cbind(nc_location, nc_adj_matrix)
for(i in 4:9) {
  nc_plot_data[, i] <- factor(nc_plot_data[, i], levels = c(">", "=", "<"))
}


{
p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = ELAVL3_vs_FAIM2), size = 0.1) +
    scale_colour_manual(values = c("red3", "gray", "blue3")) +
    theme_classic() +
    guides(color = "none") +
    ggtitle("AD", subtitle = "ELAVL3_vs_FAIM2")
p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = ELAVL3_vs_FAIM2), size = 0.1) +
  theme_classic() +
  guides(color = "none") +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  ggtitle("NC", subtitle = "ELAVL3_vs_FAIM2")
p <- p1 + p2

ggsave(filename = "results/sc_bulk/ELAVL3_vs_FAIM2.png", plot = p, device = "png", dpi = 300, width = 5, height = 2.5)
}
{
  p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = ELAVL3_vs_RPH3A), size = 0.1) +
  theme_classic() +
  guides(color = "none") +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
    guides(color = "none") +
    ggtitle("AD", subtitle = "ELAVL3_vs_RPH3A")

  p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = ELAVL3_vs_RPH3A), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
    guides(color = "none") +
    ggtitle("NC", subtitle = "ELAVL3_vs_RPH3A")
  p <- p1 + p2
ggsave(filename = "results/sc_bulk/ELAVL3_vs_RPH3A.png", plot = p, device = "png", dpi = 300, width = 5, height = 2.5)
}
{
p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = GFAP_vs_SULT4A1), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  guides(color = "none") +
  ggtitle("AD", subtitle = "GFAP_vs_SULT4A1")

p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = GFAP_vs_SULT4A1), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  guides(color = "none") +
  ggtitle("NC", subtitle = "GFAP_vs_SULT4A1")
p <- p1 + p2
ggsave(filename = "results/sc_bulk/GFAP_vs_SULT4A1.png", plot = p, device = "png", dpi = 300, width = 5, height = 2.5)
}
{
p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = LRRTM2_vs_HMP19), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  guides(color = "none") +
  ggtitle("AD", subtitle = "LRRTM2_vs_HMP19")

p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = LRRTM2_vs_HMP19), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  guides(color = "none") +
  ggtitle("NC", subtitle = "LRRTM2_vs_HMP19")
p <- p1 + p2
ggsave(filename = "results/sc_bulk/LRRTM2_vs_HMP19.png", plot = p, device = "png", dpi = 300, width = 5, height = 2.5)
}
{
p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = SLC24A2_vs_NEUROD6), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  guides(color = "none") +
  ggtitle("AD", subtitle = "SLC24A2_vs_NEUROD6")

p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = SLC24A2_vs_NEUROD6), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  guides(color = "none") +
  ggtitle("NC", subtitle = "SLC24A2_vs_NEUROD6")
p <- p1 + p2
ggsave(filename = "results/sc_bulk/NEUROD6_vs_SLC24A2.png", plot = p, device = "png", dpi = 300, width = 5, height = 2.5)
}
{
p1 <- ggplot(ad_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = RTN1_vs_SNAP25), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  guides(color = "none") +
  ggtitle("AD", subtitle = "RTN1_vs_SNAP25")
p2 <- ggplot(nc_plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = RTN1_vs_SNAP25), size = 0.1) +
  theme_classic() +
  scale_colour_manual(values = c("red3", "gray", "blue3")) +
  guides(color = "none") +
  ggtitle("NC", subtitle = "RTN1_vs_SNAP25")
p <- p1 + p2
ggsave(filename = "results/sc_bulk/RTN1_vs_SNAP25.png", plot = p, device = "png", dpi = 300, width = 5, height = 2.5)

}

{
  p1 <- ggplot(ad_plot_data,aes(x = ., y= RTN1_vs_SNAP25, fill = RTN1_vs_SNAP25)) +
    geom_bar(stat='identity') +
    guides(fill = "none") +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank())
  p2 <- ggplot(nc_plot_data,aes(x = ., y= RTN1_vs_SNAP25, fill = RTN1_vs_SNAP25)) +
    geom_bar(stat='identity') +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank()) 
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/RTN1_vs_SNAP25_bar.png", plot = p, device = "png", dpi = 300, width = 8, height = 2.5)
  
}
{
  p1 <- ggplot(ad_plot_data,aes(x = ., y= SLC24A2_vs_NEUROD6, fill = SLC24A2_vs_NEUROD6)) +
    geom_bar(stat='identity') +
    guides(fill = "none") +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank())
  p2 <- ggplot(nc_plot_data,aes(x = ., y= SLC24A2_vs_NEUROD6, fill = SLC24A2_vs_NEUROD6)) +
    geom_bar(stat='identity') +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank()) 
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/SLC24A2_vs_NEUROD6_bar.png", plot = p, device = "png", dpi = 300, width = 8, height = 2.5)
  
}
{
  p1 <- ggplot(ad_plot_data,aes(x = ., y= LRRTM2_vs_HMP19, fill = LRRTM2_vs_HMP19)) +
    geom_bar(stat='identity') +
    guides(fill = "none") +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank())
  p2 <- ggplot(nc_plot_data,aes(x = ., y= LRRTM2_vs_HMP19, fill = LRRTM2_vs_HMP19)) +
    geom_bar(stat='identity') +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank()) 
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/LRRTM2_vs_HMP19_bar.png", plot = p, device = "png", dpi = 300, width = 8, height = 2.5)
  
}
{
  p1 <- ggplot(ad_plot_data,aes(x = ., y= GFAP_vs_SULT4A1, fill = GFAP_vs_SULT4A1)) +
    geom_bar(stat='identity') +
    guides(fill = "none") +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank())
  p2 <- ggplot(nc_plot_data,aes(x = ., y= GFAP_vs_SULT4A1, fill = GFAP_vs_SULT4A1)) +
    geom_bar(stat='identity') +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank()) 
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/GFAP_vs_SULT4A1_bar.png", plot = p, device = "png", dpi = 300, width = 8, height = 2.5)
  
}
{
  p1 <- ggplot(ad_plot_data,aes(x = ., y= ELAVL3_vs_RPH3A, fill = ELAVL3_vs_RPH3A)) +
    geom_bar(stat='identity') +
    guides(fill = "none") +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank())
  p2 <- ggplot(nc_plot_data,aes(x = ., y= ELAVL3_vs_RPH3A, fill = ELAVL3_vs_RPH3A)) +
    geom_bar(stat='identity') +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank()) 
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/ELAVL3_vs_RPH3A_bar.png", plot = p, device = "png", dpi = 300, width = 8, height = 2.5)
  
}
{
  p1 <- ggplot(ad_plot_data,aes(x = ., y= ELAVL3_vs_FAIM2, fill = ELAVL3_vs_FAIM2)) +
    geom_bar(stat='identity') +
    guides(fill = "none") +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank())
  p2 <- ggplot(nc_plot_data,aes(x = ., y= ELAVL3_vs_FAIM2, fill = ELAVL3_vs_FAIM2)) +
    geom_bar(stat='identity') +
    theme_classic() +
    xlab("cell type") +
    theme(axis.text.y = element_blank()) 
  p <- p1 + p2
  ggsave(filename = "results/sc_bulk/ELAVL3_vs_FAIM2bar.png", plot = p, device = "png", dpi = 300, width = 8, height = 2.5)
  
}

