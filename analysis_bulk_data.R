# load packages -----------------------------------------------------------

source("final-code/analysis_data_functions.R", encoding = "UTF-8")

# all data co-gene ------------------------------------------------------

data_dirs <- list("data/GSE48350/GSE48350_series_matrix_ad.csv",
                  "data/GSE5281/GSE5281_series_matrix_ad.csv",
                  "data/GSE33000/GSE33000_series_matrix_ad.csv",
                  "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                  "data/GSE104704/GSE104704_RNAseq_ad.csv")

genes <- lapply(data_dirs, function(x) {
  read.csv(x, row.names = 1) %>% rownames()
  })

b_gene <- read.table("data/genes/brain_tiger.txt", sep = "\t", header = T, quote = "", fill = T) %$% 
  Gene_Symbol %>% 
  unique()
data_co_gene <- Reduce(intersect, genes) %>% intersect(b_gene)

# training set ------------------------------------------------------------

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

# training procession -----------------------------------------------------

pair_data <- xxs_get_pair_data_from_matrix(matrix_train_ad, matrix_train_nc, rate_cutoff = 0.6)
lasso_pair <- xxs_get_pair_from_lasso(pair_data, matrix_train_ad, matrix_train_nc)

lasso_pair_plot <- lasso_pair
colnames(lasso_pair_plot) <- c("Gene Pair", "Frequency")

lasso_pair_slc <- subset(lasso_pair, subset = lasso_pair$Freq > 60, select = .) %>% 
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
            "results/bulk/Active.Coefficients_train_from_bulk.txt",
            row.names = F)
{
  Coefficients <- coef(cv.fit, s = Lambda)
  Coefficients <- as.matrix(Coefficients)[-1]
  Coefficients <- data.frame(Coefficients)
  Coefficients$`Gene pair` <- rownames(Coefficients)
  ggplot(Coefficients) +
    geom_col(aes(x = `Gene pair`, y = Coefficients)) +
  theme_classic() +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) # 去掉标签
    
  xlab(NULL)# 去掉坐标名字
  lasso_pair_plot_1 <- lasso_pair_plot[order(lasso_pair_plot$Frequency, decreasing = F), ]
  lasso_pair_plot_1 <- lasso_pair_plot_1[c(150, 170, 185, 200, 215, 220, 230, 231, 234), ]
  lasso_pair_plot_1$`Gene Pair` <- factor(lasso_pair_plot_1$`Gene Pair`, levels = lasso_pair_plot_1$`Gene Pair`)
  ggplot(lasso_pair_plot_1) +
    geom_col(aes(y = `Gene Pair`, x = Frequency, fill = Frequency)) +
    # scale_fill_gradient(low = "black", high = "#FFA54F") +
    # geom_hline(yintercept = -60, col = "red", linetype = 2) +
    scale_fill_gradient(low = "lightpink", high= "red3") +
    
    # xlim(30, 100) +
    theme_classic() +
    # theme(axis.line = element_line(colour = "white")) + # 坐标轴白色
    # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # 去掉标签
    # theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + # 去掉标签
    # xlab(NULL)# 去掉坐标名字
    ylab(NULL) # 去掉坐标名字
  
  }
{
  ggplot(lasso_pair_plot) +
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
  
  ggsave("results/bulk/gene_pair_freq.png", device = "png", width = 8, height = 8, dpi = 300)
} # plot_code
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
  ggsave("results/bulk/slc_pair_coeff.png", device = "png", width = 8, height = 8, dpi = 300)
  ggsave("results/bulk/slc_pair_coeff_1.png", device = "png", width = 5, height = 3, dpi = 300)
}
# plot_code

pred_train = predict(cv.fit, matrix_cbind_train, type = "response", s = Lambda) %>%  #s=seLambda
  as.numeric()
{
  png("results/bulk/bulk pair data train auc.png", width = 900, height = 900, units = "px", res = 300)
auc1 <- roc(label_train, pred_train, main = "Training",
            levels = c("0", "1"), plot = T, print.auc = T, direction = "<", legacy.axes = T)
dev.off()
}
# plot(auc1, main = round(auc1$auc, 4))


pred_test = predict(cv.fit, matrix_cbind_test, type="response", s = Lambda) %>%  #s=seLambda
  as.numeric()
{
png("results/bulk/bulk pair data test auc.png", width = 900, height = 900, units = "px", res = 300)
auc2 <- roc(label_test, pred_test, main = "Testing",
            levels = c("0", "1"), plot = T, print.auc = T, direction = "<", legacy.axes = T)
dev.off()
}
# plot(auc2, legacy.axes = T, main = round(auc2$auc, 4))

ad_gene <- rownames(Active.Coefficients) %>%
  str_split(., "_vs_", simplify = T) %>%
  c() %>% 
  unique()
write.table(ad_gene, 
            "results/bulk/gene_train_from_bulk.txt",
            row.names = F)
ad_gene <- read.table("results/bulk/gene_train_from_bulk.txt",
                      header = T) %>%
  .[, 1]
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit,
                        "results/bulk/", "GSE48350", 
                        "data/GSE48350/GSE48350_series_matrix_ad.csv",
                        "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit,
                        "results/bulk/", "GSE5281",
                        "data/GSE5281/GSE5281_series_matrix_ad.csv",
                        "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit,
                        "results/bulk/", "GSE33000",
                        "data/GSE33000/GSE33000_series_matrix_ad.csv",
                        "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit,
                        "results/bulk/", "GSE104704", 
                        "data/GSE104704/GSE104704_RNAseq_ad.csv",
                        "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_auc_from_lasso_pair(lasso_pair_slc, cv.fit,
                        "results/bulk/", "GSE159699", 
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

RP_inter <- read.table("../iscience/Download_data_RP.txt", header = T, sep = "\t", fill = T, quote = "")
RP_inter <- subset(RP_inter,subset =
                     RP_inter$Species1 == "Homo sapiens"
                   & RP_inter$Species2 == "Homo sapiens"
                   & RP_inter$Interactor1.Symbol %in% ad_gene
                   & RP_inter$Interactor2.Symbol %in% ad_gene)
nrow(RP_inter)

RP_inter <- RP_inter[, c(2, 5)]
RP_inter$Interactor2.Symbol <- paste0(RP_inter$Interactor2.Symbol, " ")

links <- RP_inter

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
  png("results/bulk/net.png", width = 1500, height = 1500, res = 300)
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
write.csv(netdata, "results/bulk/netdata.csv", row.names = F)

netdata$ad_gene <- factor(netdata$ad_gene, levels = netdata$ad_gene)
{
  ggplot(netdata, aes(y = ad_gene, x = links, fill = -links)) +
    geom_col() +
    ylab(label = "Gene") +
    xlab(label = "The number of edges") +
    guides(fill = "none") +
    scale_fill_distiller(palette = "Reds") +
    theme_classic()
ggsave("results/bulk/gene_net_freq.png", width = 3, height = 3, device = "png")
}
# 染色体 ---------------------------------------------------------------------
{

png("results/bulk/chome_gene.png", width = 1200, height = 1200, res = 300)
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
track.num <- 1

CTNND2 <- c("chr5", 10971836, 11904446, "CTNND2")
ELAVL3 <- c("chr19", 11451326, 11481046, "ELAVL3")
KCNQ2 <- c("chr20", 63400208, 63472655, "KCNQ2")
FAM19A5 <- c("chr1", 17271103, 17692943, "FAM19A5")
PSD2 <- c("chr5", 139742475, 139844466, "PSD2")
SLC39A12 <- c("chr10", 17951918, 18043285, "SLC39A12")
TNR <- c("chr1", 175315194, 175743595, "TNR")
ASB16 <- c("chr17", 44170704, 44179084, "ASB16")
C1orf61 <- c("chr1", 156404252, 156429548, "C1orf61")
MOG <- c("chr6", 29657092, 29672365, "MOG")
RCircos.Gene.Label.Data_1 <- data.frame(CTNND2, 
                                        ELAVL3,
                                        KCNQ2,
                                        FAM19A5,
                                        PSD2,
                                        SLC39A12,
                                        TNR,
                                        ASB16,
                                        C1orf61,
                                        MOG) %>%
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

RCircos.Link.Data_1 <- data.frame(c(TNR, ASB16),
                                  c(SLC39A12, MOG),
                                  c(PSD2, C1orf61),
                                  c(PSD2, ASB16),
                                  c(FAM19A5, KCNQ2),
                                  c(KCNQ2, ELAVL3),
                                  c(CTNND2, ASB16)) %>%
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
                       pvalueCutoff = 0.1,
                       qvalueCutoff = 0.1)
all_enrich <- setReadable(all_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
png("results/bulk/go-enrich-net.png", width = 2700, height = 2100, res = 300)

enrichplot::cnetplot(all_enrich,
                     circular = FALSE,
                     showCategory = 10,
                     layout = "circle",
                     # node_label="gene",
                     colorEdge = TRUE
                     )
dev.off()

all_enrich <- enrichGO(geneid$ENTREZID,
                       OrgDb = org.Hs.eg.db, 
                       ont = "ALL",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1)

png("results/bulk/all.png", width = 2700, height = 3600, res = 300)
barplot(all_enrich, drop = TRUE, showCategory = 10, split = "ONTOLOGY", ) + 
  facet_grid(ONTOLOGY~., scale = 'free')
dev.off()

bp_enrich <- enrichGO(geneid$ENTREZID,
                      OrgDb = org.Hs.eg.db, 
                       ont = "BP",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1)
png("results/bulk/bp.png", width = 1800, height = 1800, res = 300)
dotplot(bp_enrich, title = "BP")
dev.off()
cc_enrich <- enrichGO(geneid$ENTREZID,
                      OrgDb = org.Hs.eg.db, 
                       ont = "CC",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1)
png("results/bulk/cc.png", width = 1800, height = 1800, res = 300)
dotplot(cc_enrich, title = "CC")
dev.off()
mf_enrich <- enrichGO(geneid$ENTREZID,
                      OrgDb = org.Hs.eg.db, 
                       ont = "MF",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1)
png("results/bulk/mf.png", width = 1800, height = 1800, res = 300)
dotplot(mf_enrich, title = "MF")
dev.off()

# both of enrich have no results
# don't need to write functions

# expression level --------------------------------------------------------

# percent bar plot
xxs_exp_col_plot(lasso_pair_slc, "results/bulk/", "GSE48350", 
                 "data/GSE48350/GSE48350_series_matrix_ad.csv",
                 "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_exp_col_plot(lasso_pair_slc, "results/bulk/", "GSE5281",
                 "data/GSE5281/GSE5281_series_matrix_ad.csv",
                 "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_exp_col_plot(lasso_pair_slc, "results/bulk/", "GSE33000",
                 "data/GSE33000/GSE33000_series_matrix_ad.csv",
                 "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_exp_col_plot(lasso_pair_slc, "results/bulk/", "GSE104704", 
                 "data/GSE104704/GSE104704_RNAseq_ad.csv",
                 "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_exp_col_plot(lasso_pair_slc, "results/bulk/", "GSE159699", 
                 "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                 "data/GSE159699/GSE159699_RNAseq_count_nc.csv")


# waterfall --------------------------------------------------------------

xxs_waterfall_7(lasso_pair_slc, "GSE48350", "results/bulk/", 
                "data/GSE48350/GSE48350_series_matrix_ad.csv",
                "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_waterfall_7(lasso_pair_slc, "GSE5281", "results/bulk/", 
                "data/GSE5281/GSE5281_series_matrix_ad.csv",
                "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_waterfall_7(lasso_pair_slc, "GSE33000", "results/bulk/", 
                "data/GSE33000/GSE33000_series_matrix_ad.csv",
                "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_waterfall_7(lasso_pair_slc, "GSE104704", "results/bulk/", 
                "data/GSE104704/GSE104704_RNAseq_ad.csv",
                "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_waterfall_7(lasso_pair_slc, "GSE159699", "results/bulk/",  
                "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
                "data/GSE159699/GSE159699_RNAseq_count_nc.csv")

# boxplot ---------------------------------------------------------------

xxs_boxplot(lasso_pair_slc, "results/bulk/", "GSE48350", 
            "data/GSE48350/GSE48350_series_matrix_ad.csv",
            "data/GSE48350/GSE48350_series_matrix_nc.csv")
xxs_boxplot(lasso_pair_slc, "results/bulk/", "GSE5281",
            "data/GSE5281/GSE5281_series_matrix_ad.csv",
            "data/GSE5281/GSE5281_series_matrix_nc.csv")
xxs_boxplot(lasso_pair_slc, "results/bulk/", "GSE33000",
            "data/GSE33000/GSE33000_series_matrix_ad.csv",
            "data/GSE33000/GSE33000_series_matrix_nc.csv")
xxs_boxplot(lasso_pair_slc, "results/bulk/", "GSE104704", 
            "data/GSE104704/GSE104704_RNAseq_ad.csv",
            "data/GSE104704/GSE104704_RNAseq_nc.csv")
xxs_boxplot(lasso_pair_slc, "results/bulk/", "GSE159699", 
            "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
            "data/GSE159699/GSE159699_RNAseq_count_nc.csv")