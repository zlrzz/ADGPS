# ADGPS

**Alzheimer's Disease Gene Pair Signature: single-cell-informed classification of bulk brain transcriptomes.**

ADGPS represents expression profiles as within-sample gene comparisons and uses LASSO logistic regression to distinguish Alzheimer's disease (AD) from non-AD controls (NC). This repository contains R scripts for data preparation, a bulk-only baseline, single-cell feature selection, transfer of selected pairs to bulk data, and visualization.


## Requirements and installation

Use R with a consistent CRAN/Bioconductor environment. RStudio is optional; opening `AD.Rproj` is convenient for interactive execution. Original package versions are not pinned in `final-code/`, so check compatibility before a full run and record `sessionInfo()`.

Install the packages imported by the scripts and supporting dependencies in R:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
cran_packages <- c(
  "glmnet", "RCircos", "magrittr", "stringr", "ggplot2", "pROC",
  "igraph", "RColorBrewer", "MetBrewer", "VennDiagram", "ggalluvial",
  "pheatmap", "ggExtra", "Seurat", "downloader", "rvest", "hdf5r",
  "BiocManager"
)
install.packages(cran_packages)
BiocManager::install(
  c("clusterProfiler", "org.Hs.eg.db", "glmGamPoi", "GEOquery", "enrichplot"),
  ask = FALSE, update = FALSE
)
```

The original installer omits explicit installation of 'downloader', 'rvest', 'GEOquery', 'enrichplot', and the HDF5 reader dependency. See the [Bioconductor installation guide](https://bioconductor.org/install/) for R/package compatibility. Single-cell scripts use Seurat's [`Read10X_h5()`](https://satijalab.org/seurat/reference/read10x_h5) to read Cell Ranger HDF5 matrices.

The shared analysis helper imports the modeling and visualization packages when sourced. Even a bulk-only run currently loads Seurat and other visualization dependencies.

## Data preparation

### Datasets used by the scripts

| GEO accession | Role in this implementation | Processed input |
| --- | --- | --- |
| [GSE129308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129308) | Single-cell feature selection and Seurat analysis | Named `*_filtered_feature_bc_matrix.h5` files |
| [GSE48350](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350) | Pooled bulk development and cohort evaluation | `GSE48350_series_matrix_{ad,nc}.csv` |
| [GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281) | Pooled bulk development and cohort evaluation | `GSE5281_series_matrix_{ad,nc}.csv` |
| [GSE33000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33000) | Additional bulk cohort evaluation | `GSE33000_series_matrix_{ad,nc}.csv` |
| [GSE159699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159699) | Additional bulk cohort evaluation | `GSE159699_RNAseq_count_{ad,nc}.csv` |
| [GSE104704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104704) | Additional bulk cohort evaluation | `GSE104704_RNAseq_{ad,nc}.csv` |

`{ad,nc}` denotes two separate filenames. Roles reflect the code's usage. GSE48350 and GSE5281 are reused for cohort evaluation after pooling; those evaluations are not independent external validation.

Bulk CSVs contain gene symbols in the first column and numeric expression values in sample columns, with one row per gene. They are read using `read.csv(..., row.names = 1)`. Before combining matrices, check unique gene names, matching row order, and AD/NC group assignments.

### Download and preprocessing

Run from the repository root, which contains `final-code/`, `data/`, and `results/`:

```r
stopifnot(dir.exists("final-code"))
dir.create("data", showWarnings = FALSE)
for (path in c("results/bulk", "results/sc", "results/sc_bulk")) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
source("final-code/dowmload_data.R", encoding = "UTF-8")
```

Supply the following before preprocessing or single-cell analysis:

- `data/GSE*/sample.csv` metadata exports where the script reads them. Fields include `Accession`, `Title`, and, for the GSE33000 grouping block, `Sample.Type`. Check each dataset's grouping code against its metadata.
- `GPL570-55999.txt` under both `data/GSE5281/` and `data/GSE48350/`, with a `Gene.Symbol` annotation column.
- `data/genes/brain_tiger.txt`, a tab-delimited table containing `Gene_Symbol`.
- HDF5 matrices under `data/GSE129308/`. The combined workflow uses `GSM3704357_1-MAP2_filtered_feature_bc_matrix.h5` and `GSM6261344_Control-1-MAP2_filtered_feature_bc_matrix.h5`. Other variants reference additional sample files.

The downloader fetches bulk expression files and prints metadata reminders. It does not download these additional inputs or verify their completeness.

```r
source("final-code/process_data.R", encoding = "UTF-8")
```

Preprocessing removes blank gene symbols and averages expression rows sharing a gene symbol. Groups are selected through dataset-specific metadata or sample-name patterns. `process_data.R` clears the global environment at completion, so initialize analysis objects afterward.

## Running the analysis

### 1. Load helpers and define the common genes

This reproduces the initialization in `analysis_bulk_data.R` without fitting the bulk model:

```r
source("final-code/analysis_data_functions.R", encoding = "UTF-8")
data_dirs <- c(
  "data/GSE48350/GSE48350_series_matrix_ad.csv",
  "data/GSE5281/GSE5281_series_matrix_ad.csv",
  "data/GSE33000/GSE33000_series_matrix_ad.csv",
  "data/GSE159699/GSE159699_RNAseq_count_ad.csv",
  "data/GSE104704/GSE104704_RNAseq_ad.csv"
)
stopifnot(all(file.exists(data_dirs)))
genes <- lapply(data_dirs, function(path) rownames(read.csv(path, row.names = 1)))
brain <- read.table(
  "data/genes/brain_tiger.txt", header = TRUE, sep = "\t",
  quote = "", fill = TRUE
)
b_gene <- unique(brain$Gene_Symbol)
data_co_gene <- intersect(Reduce(intersect, genes), b_gene)
stopifnot(length(data_co_gene) > 1L)
```

The source defines common genes using gene availability across the five listed bulk AD matrices and the brain gene list. The universe is not derived solely from the pooled training partition.

### 2. Fit the bulk-only baseline

Open `final-code/analysis_bulk_data.R` and execute the model-fitting and validation sections, stopping before `# net` unless the additional network resources are available. It pools GSE48350 and GSE5281, uses seed 38 for separate 70/30 AD and NC sample splits, screens pairs, repeats LASSO selection, fits the classifier, and writes coefficient tables and ROC plots under `results/bulk/`.

### 3. Fit the single-cell model and transfer features to bulk

After initializing `data_co_gene`, open `final-code/addition_sce.R` and run its sections in order:

1. Load the named AD/NC matrices and perform SCTransform integration, PCA, UMAP, and clustering.
2. Inspect marker plots and review cluster-to-cell-type assignments. Cluster identifiers are hard-coded.
3. Run single-cell pair screening, repeated selection, fitting, and internal evaluation.
4. Run `transfer learning` and validation to refit the selected pairs on pooled bulk data and write outputs under `results/sc_bulk/`.
5. Execute downstream plotting when its additional resources are available.

The combined script splits cells from the selected AD/NC files 70/30. This is cell-level internal evaluation, not held-out-donor evaluation. Alternative scripts use different sample selections and splits and should be treated as separate workflows.

### 4. Additional analyses

Network sections read `../iscience/Download_data_RR.txt` and `../iscience/Download_data_RP.txt`, which the downloader does not supply. Prepare these tables or update their paths. Chromosome diagrams, six-/seven-pair waterfall plots, and some single-cell overlays use fixed gene lists or model-specific assumptions that need review after refitting.

`final-code/test.R` compares bulk-only and single-cell-informed models across 100 outer random splits and writes `results/auc_all.csv`. It requires objects from earlier sections, and its ROC helper calls need alignment with the current helper interface.

## Outputs

| Location | Contents |
| --- | --- |
| `data/GSE*/` | Processed AD/NC CSVs and preprocessing logs (`output.txt`). |
| `results/bulk/` | Bulk-only coefficient/gene tables, internal/cohort ROC plots, and downstream figures. |
| `results/sc/` | Single-cell model plots, integration figures, and cell-type visualizations, depending on the variant. |
| `results/sc_bulk/` | Transferred model coefficient/gene tables, bulk ROC plots, enrichment and expression figures. |
| `results/auc_all.csv` | Repeated comparison results from `test.R` when its prerequisites are satisfied. |

Bulk and transfer scripts write `Active.Coefficients_train_from_bulk.txt` and `gene_train_from_bulk.txt` in their respective directories. The coefficient tables omit the intercept and are not complete serialized prediction models. Save the fitted object, selected pair order, and chosen lambda to reproduce predictions.

Existing outputs may come from different variants or earlier runs. Refitting can change the selected pairs; plotting helpers do not guarantee that a new fit selects exactly six or seven pairs.


## Reproducibility

Record inputs, grouping decisions, common genes, selected pair order, seeds, package versions, and the lambda used in each evaluation. Preserve previous results before rerunning because scripts write to fixed filenames.

For example, after fitting the transferred bulk model in the combined workflow:

```r
saveRDS(
  list(
    fit = cv.fit,
    pairs = lasso_pair_slc_sce_to_bulk,
    lambda = Lambda,
    common_genes = data_co_gene,
    session_info = sessionInfo()
  ),
  file = "results/sc_bulk/ADGPS_model.rds"
)
```

This README was checked against the current source. Syntax inspection is separate from a complete analysis run; dependency installation, all-input download, and reproduction of every figure were not performed when preparing this documentation.
