# functions for process data 
# xxs_rownames_by_gene ----------------------------------------------------

xxs_rownames_by_gene <- function(gene_symbol, data) {
  
  if("" %in% gene_symbol) {
    index_del <- which(gene_symbol == "")
    gene_symbol <- gene_symbol[-index_del]
    data <- data[-index_del, ]
    print(paste0(length(index_del), " blank gene symbol had beed deleted"))
  }
  
  if(T %in% duplicated(gene_symbol)) {
    gene_symbol <- factor(gene_symbol)
    data <- apply(data, 2, function(x) {tapply(x, gene_symbol, mean)}) %>% # average the expression by the same gene
      data.frame()
    print("Average duplicate genes")
  } else {
    rownames(data) <- gene_symbol
  }

  return(data)
}

# xxs_sample_group --------------------------------------------------------

xxs_sample_group <- function(sample, chr_1, chr_2) {
  
  library(magrittr) # 管道编辑工具
  library(stringr) # 字符编辑
  library(rvest)  # 网页抓取工具
  
  for (i in 1:nrow(sample)) {
    url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", sample[i, 1])
    HTML <- read_html(x = url)
    sample[i, 3] <- html_text(HTML) %>% 
      str_split(., chr_1, simplify = T) %>%
      .[2] %>%
      str_split(., chr_2, simplify = T) %>%
      .[1]
  }
  
  return(sample)
}

# End ---------------------------------------------------------------------

print("functions for process_data load completed")
