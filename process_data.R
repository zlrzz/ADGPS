# packages ----------------------------------------------------------------

library(stringr)
library(magrittr)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
source("final-code/process_data_functions.R", encoding = "UTF-8")

# GSE5281 -----------------------------------------------------------------

sink("data/GSE5281/output.txt", split = T)

gene_symbol <- read.table("data/GSE5281/GPL570-55999.txt", 
                          header = T, sep = "\t", fill = T, quote = "") %$% 
  Gene.Symbol %>% 
  str_split(., "///", simplify = T) %>% 
  .[, 1]

series_matrix <- read.table("data/GSE5281/GSE5281_series_matrix.txt.gz",
                            header = T, sep = "\t", fill = T, quote = "", comment.char = "!", check.names = F, row.names = 1) %>% 
  xxs_rownames_by_gene(gene_symbol, .)

sample <- read.csv("data/GSE5281/sample.csv")

temp_function <- function(chr, dir) {
  grep(chr, sample$Title) %>% 
    sample$Accession[.] %>% 
    paste0("X.", ., ".") %>% 
    match(., colnames(series_matrix)) %>% 
    series_matrix[, .] %T>% 
    write.csv(., dir) %>% 
    ncol() %>% 
    paste0(., " - ", chr) %>% 
    print()
}

temp_function("affected", "data/GSE5281/GSE5281_series_matrix_ad.csv")
temp_function("control", "data/GSE5281/GSE5281_series_matrix_nc.csv")

print("files in the dir are: ")
print(list.files("data/GSE5281/"))
sink()

# GSE48350 ----------------------------------------------------------------

sink("data/GSE48350/output.txt", split = T)

gene_symbol <- read.table("data/GSE48350/GPL570-55999.txt", 
                          header = T, sep = "\t", fill = T, quote = "") %$% 
  Gene.Symbol %>% 
  str_split(., "///", simplify = T) %>% 
  .[, 1]

series_matrix <- read.table("data/GSE48350/GSE48350_series_matrix.txt.gz",
                            header = T, sep = "\t", fill = T, quote = "", comment.char = "!", check.names = F, row.names = 1) %>% 
  xxs_rownames_by_gene(gene_symbol, .)

sample <- read.csv("data/GSE48350/sample.csv")

temp_function <- function(chr, dir) {
  grep(chr, sample$Title) %>% 
    sample$Accession[.] %>% 
    paste0("X.", ., ".") %>% 
    match(., colnames(series_matrix)) %>% 
    series_matrix[, .] %T>% 
    write.csv(., dir) %>% 
    ncol() %>% 
    paste0(., " - ", chr) %>% 
    print()
}

temp_function("AD", "data/GSE48350/GSE48350_series_matrix_ad.csv")
temp_function("_indiv", "data/GSE48350/GSE48350_series_matrix_nc.csv")

print("files in the dir are: ")
print(list.files("data/GSE48350/"))
sink()

# GSE33000 ----------------------------------------------------------------

sink("data/GSE33000/output.txt", split = T)

data <- read.table("data/GSE33000/GSE33000_raw_data.txt.gz",
                   header = T, sep = "\t", fill = T, quote = "", check.names = F)
gene_symbol <- data$Gene %>% 
  str_split(., ",", simplify = T) %>% 
  .[, 1] %>% 
  gsub('\"', "", .)

sample <- read.csv("data/GSE33000/sample.csv") %>% 
  xxs_sample_group(., "disease status: ", "\nTreatment protocol")

series_matrix <- read.table("data/GSE33000/GSE33000_series_matrix.txt.gz",
                            header = T, sep = "\t", fill = T, quote = "", comment.char = "!", check.names = F, row.names = 1) %>% 
  xxs_rownames_by_gene(gene_symbol, .)

temp_function <- function(chr, dir) {
  grep(chr, sample$Sample.Type) %>% 
    sample$Accession[.] %>% 
    paste0("X.", ., ".") %>% 
    match(., colnames(series_matrix)) %>% 
    series_matrix[, .] %T>% 
    write.csv(., dir) %>% 
    ncol() %>% 
    paste0(., " - ", chr) %>% 
    print()
}

temp_function("Alzheimer's disease", "data/GSE33000/GSE33000_series_matrix_ad.csv")
temp_function("non-demented", "data/GSE33000/GSE33000_series_matrix_nc.csv")

print("files in the dir are: ")
print(list.files("data/GSE33000/"))
sink()

# GSE159699 ---------------------------------------------------------------

sink("data/GSE159699/output.txt", split = T)
data <- read.table("data/GSE159699/GSE159699_summary_count.star.txt.gz",
                   header = T, sep = "\t", fill = T, quote = "")
gene_symbol <- data$refGene
data <- data[, -1] %>% 
  xxs_rownames_by_gene(gene_symbol, .)

temp_function <- function(chr, dir) {
  grep(chr, colnames(data)) %>% 
    data[, .] %T>% 
    write.csv(., dir) %>% 
    ncol() %>% 
    paste0(., " - ", chr) %>% 
    print()
}

temp_function("AD", "data/GSE159699/GSE159699_RNAseq_count_ad.csv")
temp_function("Old", "data/GSE159699/GSE159699_RNAseq_count_nc.csv")

print("files in the dir are: ")
print(list.files("data/GSE159699/"))
sink()

# GSE104704 ---------------------------------------------------------------

sink("data/GSE104704/output.txt", split = T)
data <- read.table("data/GSE104704/GSE104704_RNA-Seq_Table.txt.gz",
                   header = T, sep = "\t", fill = T, quote = "", check.names = F)
gene_symbol <- data$Gene
data <- data[, -c(1, 32, 33, 34)] %>% 
  xxs_rownames_by_gene(gene_symbol, .)

sample <- read.csv("data/GSE104704/sample.csv")
temp_function <- function(chr, dir) {
  match(str_split(colnames(data), ".RNA", simplify = T)[, 1], 
        str_split(sample$Title, ".RNA", simplify = T)[, 1]) %>% 
    sample$Title[.] %>%
    grep(chr, .) %>% 
    data[, .] %T>% 
    write.csv(., dir) %>% 
    ncol() %>% 
    paste0(., " - ", chr) %>% 
    print()
}

temp_function("AD", "data/GSE104704/GSE104704_RNAseq_ad.csv")
temp_function("Old", "data/GSE104704/GSE104704_RNAseq_nc.csv")

print("files in the dir are: ")
print(list.files("data/GSE104704/"))
sink()

# end ---------------------------------------------------------------------

rm(list = ls())
gc()
