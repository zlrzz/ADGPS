# load packages -----------------------------------------------------------

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
sce_ad_train <- sce_ad[, train_index_ad]
sce_nc_train <- sce_nc[, train_index_nc]

# 内部拆分为训练和验证组 -------------------------------------------------------------

pair_data_sce <- xxs_get_pair_data_from_matrix(sce_ad_train, sce_nc_train, rate_cutoff = 0.6)
lasso_pair_sce <- xxs_get_pair_from_lasso(pair_data_sce, sce_ad_train, sce_nc_train)

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

ggplot(Active.Coefficients_sce[1:10, ]) +
  geom_point(aes(x = Active.Coefficients_sce, y = `Gene Pair`)) +
  theme(axis.title.y = NA) +
  theme_classic()

pred_train <- predict(sce_cv.fit, sce_matrix_cbind_train, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
auc1 <- roc(sce_label_train, pred_train, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Training")


# sce_ad_test <- Read10X_h5("data/GSE129308/GSM3704371_8-MAP2_filtered_feature_bc_matrix.h5") %>%
#   data.frame() %>%
#   xxs_filter_matrix_by_gene(gene = data_co_gene)
# sce_nc_test <- Read10X_h5("data/GSE129308/GSM6261351_Control-8-MAP2_filtered_feature_bc_matrix.h5") %>%
#   data.frame() %>%
#   xxs_filter_matrix_by_gene(gene = data_co_gene)

sce_label_test <- c(rep(1, ncol(sce_ad_test)), rep(0, ncol(sce_nc_test))) 
sce_matrix_cbind_test <- xxs_get_adj_matrix(lasso_pair_sce_slc, sce_ad_test, sce_nc_test)

pred_test <- predict(sce_cv.fit, sce_matrix_cbind_test, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
auc2 <- roc(sce_label_test, pred_test, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Testing-8")

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

lasso_pair_bulk <- xxs_get_pair_from_lasso(pair_data_from_sce, matrix_train_ad, matrix_train_nc)

lasso_pair_plot <- subset(lasso_pair_bulk, subset = lasso_pair_bulk$Freq > 30)
lasso_pair_plot <- lasso_pair_plot[order(lasso_pair_plot$Freq, decreasing = F), ]
lasso_pair_plot$. <- factor(lasso_pair_plot$., levels = lasso_pair_plot$.)
colnames(lasso_pair_plot) <- c("Gene Pair", "Freq")
ggplot(lasso_pair_plot) +
  geom_col(aes(x = Freq, y = `Gene Pair`)) +
  theme_classic()

lasso_pair_slc <- subset(lasso_pair_bulk, subset = lasso_pair_bulk$Freq > 60, select = .) %>% 
  .[, 1] %>% 
  str_split(., "_vs_", simplify = T)

label_train <- c(rep(1, ncol(matrix_train_ad)), rep(0, ncol(matrix_train_nc))) 
matrix_cbind_train <- xxs_get_adj_matrix(lasso_pair_slc, matrix_train_ad, matrix_train_nc)

label_test <- c(rep(1, ncol(matrix_test_ad)), rep(0, ncol(matrix_test_nc))) 
matrix_cbind_test <- xxs_get_adj_matrix(lasso_pair_slc, matrix_test_ad, matrix_test_nc)

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
            "results/Active.Coefficients_train_from_bulk.txt",
            row.names = F)

ggplot(Active.Coefficients) +
  geom_point(aes(x = Active.Coefficients, y = `Gene Pair`)) +
  theme_classic()

pred_train <- predict(cv.fit, matrix_cbind_train, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
auc1 <- roc(label_train, pred_train, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Training")
# plot(auc1, main = round(auc1$auc, 4))


pred_test <- predict(cv.fit, matrix_cbind_test, type="response", s = Lambda) %>%  #s=seLambda
  as.numeric()
auc2 <- roc(label_test, pred_test, levels = c("0", "1"), plot = T, print.auc = T, direction = "<", main = "Testing")
# plot(auc2, main = round(auc2$auc, 4))

ad_gene <- rownames(Active.Coefficients) %>%
  str_split(., "_vs_", simplify = T) %>%
  c() %>% 
  unique()
write.table(ad_gene, 
            "results/gene_train_from_bulk.txt",
            row.names = F)

# validation --------------------------------------------------------------

xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit, "GSE48350", 
                        "data/GSE48350/GSE48350_series_matrix_ad.csv",
                        "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit, "GSE5281",
                        "data/GSE5281/GSE5281_series_matrix_ad.csv",
                        "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit, "GSE33000",
                        "data/GSE33000/GSE33000_series_matrix_ad.csv",
                        "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit, "GSE104704", 
                        "data/GSE104704/GSE104704_RNAseq_ad.csv",
                        "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit, "GSE159699", 
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
plot(network, 
     vertex.size = 24, 
     vertex.label.cex = 0.7, 
     vertex.label.dist = 3,
     vertex.color = my_color)

netdata <- data.frame(ad_gene)

netdata$links <- 0
for(i in 1:nrow(netdata)) {
  n <- grep(paste0("^", netdata$ad_gene[i]), c(c(links$Interactor1.Symbol), c(links$Interactor2.Symbol))) %>% 
    length()
  netdata$links[i] <- n
}
netdata <- netdata[order(netdata$links, decreasing = F), ]
write.csv(netdata, "results/netdata_sce.csv", row.names = F)

netdata$ad_gene <- factor(netdata$ad_gene, levels = netdata$ad_gene)
ggplot(netdata) +
  geom_col(aes(y = ad_gene, x = links)) +
  theme_classic()

# 染色体 ---------------------------------------------------------------------

library(RCircos)
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



CSPG5 <- c("chr3", 47562238, 47580240, "CSPG5")
ELAVL3 <- c("chr19", 11451326, 11481046, "ELAVL3")
GFAP <- c("chr17", 44903159, 44915500, "GFAP")
NEUROD6 <- c("chr7", 31337465, 31340726, "NEUROD6")
PSD2 <- c("chr5", 139742475, 139844466, "PSD2")
RTN1 <- c("chr14", 59595976, 59870776, "RTN1")
ZNF365 <- c("chr10", 62374369, 62480285, "ZNF365")
SULT4A1 <- c("chr22", 43824509, 43862513, "SULT4A1")
FAIM2 <- c("chr12", 49866896, 49903900, "FAIM2")
RPH3A <- c("chr12", 112575236, 112898881, "RPH3A")
SLC24A2 <- c("chr9", 19507455, 20307892, "SLC24A2")
HMP19 <- c("chr5", 174045706, 174109179, "HMP19")
SNAP25 <- c("chr20", 10218830, 10307418, "SNAP25")

RCircos.Gene.Label.Data_1 <- data.frame(CSPG5,
                                        ELAVL3,
                                        GFAP,
                                        NEUROD6,
                                        PSD2,
                                        RTN1,
                                        ZNF365,
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

RCircos.Link.Data_1 <- data.frame(c(ZNF365, RPH3A),
                                  c(ELAVL3, FAIM2),
                                  c(GFAP, SULT4A1),
                                  c(RTN1, SNAP25),
                                  c(ELAVL3, RPH3A),
                                  c(PSD2, HMP19),
                                  c(CSPG5, SULT4A1),
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
dotplot(bp_enrich, title = "BP")
cc_enrich <- enrichGO(geneid$ENTREZID,
                      OrgDb = org.Hs.eg.db, 
                      ont = "CC",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)
dotplot(cc_enrich, title = "CC")
mf_enrich <- enrichGO(geneid$ENTREZID,
                      OrgDb = org.Hs.eg.db, 
                      ont = "MF",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)
dotplot(mf_enrich, title = "MF")

# both of enrich have no results
# don't need to write functions

# expression level --------------------------------------------------------

lasso_pair_slc[5, ] <- lasso_pair_slc[5, 2:1]
lasso_pair_slc
# percent bar plot
xxs_exp_col_plot("GSE48350", 
                 "data/GSE48350/GSE48350_series_matrix_ad.csv",
                 "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_exp_col_plot("GSE5281",
                 "data/GSE5281/GSE5281_series_matrix_ad.csv",
                 "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_exp_col_plot("GSE33000",
                 "data/GSE33000/GSE33000_series_matrix_ad.csv",
                 "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_exp_col_plot("GSE104704", 
                 "data/GSE104704/GSE104704_RNAseq_ad.csv",
                 "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_exp_col_plot("GSE159699", 
                 "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                 "data/GSE159699/GSE159699_RNAseq_count_nc.csv")


# waterfall --------------------------------------------------------------

xxs_waterfall("GSE48350", 
              "data/GSE48350/GSE48350_series_matrix_ad.csv",
              "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_waterfall("GSE5281",
              "data/GSE5281/GSE5281_series_matrix_ad.csv",
              "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_waterfall("GSE33000",
              "data/GSE33000/GSE33000_series_matrix_ad.csv",
              "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_waterfall("GSE104704", 
              "data/GSE104704/GSE104704_RNAseq_ad.csv",
              "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_waterfall("GSE159699", 
              "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
              "data/GSE159699/GSE159699_RNAseq_count_nc.csv")

# boxplot -----------------------------------------------------------------

GSE159699_AD <- read.csv("data/GSE159699/GSE159699_RNAseq_count_ad.csv", row.names = 1) %>% xxs_filter_matrix_by_gene(., ad_gene)
GSE159699_NC <- read.csv("data/GSE159699/GSE159699_RNAseq_count_nc.csv", row.names = 1) %>% xxs_filter_matrix_by_gene(., ad_gene)

boxplot(t(GSE159699_AD[6:7, ]), main = "AD")
boxplot(t(GSE159699_NC[6:7, ]), main = "NC")
