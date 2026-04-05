README: Script 1-4 

# GSE63577 RNA‑seq analysis in R

This repository contains R scripts to download, preprocess, normalize, analyze and functionally interpret bulk RNA‑seq data from the GEO series **GSE63577: “RNA‑seq of human fibroblasts during replicative senescence”**.

The pipeline is organized into four main scripts:

- `GSE63577_preprocessing.R` – download and filter raw count data.
- `GSE63577_normalization_DEG.R` – normalization (TMM + DESeq2 vst), QC, and differential expression.
- `GSE63577_postprocessing_visualization.R` – post‑processing of DEGs, lncRNA/protein‑coding overlaps, UpSet, volcano plots, and heatmaps.
- `GSE63577_enrichment.R` – GO/KEGG functional enrichment (g:Profiler) and summary plots across cell lines.

---

## Dataset

- **Accession:** GSE63577  
- **Source:** NCBI Gene Expression Omnibus (GEO)  
- **Organism:** *Homo sapiens*  
- **Description:** RNA‑seq of human fibroblasts (BJ, WI‑38, IMR‑90, HFF, MRC‑5) at young and replicative‑senescent passages.

Raw counts are downloaded directly from GEO:

- `GSE63577_raw_counts_GRCh38.p13_NCBI.tsv.gz`

---

## Script 1 – Preprocessing and filtering (`GSE63577_preprocessing.R`)

This script:

1. **Loads dependencies**  
   Uses common Bioconductor packages for bulk RNA‑seq (e.g. `edgeR`, `DESeq2`, `SummarizedExperiment`, `ComplexHeatmap`).

2. **Downloads and imports raw counts**  
   - Builds the GEO download URL for `GSE63577_raw_counts_GRCh38.p13_NCBI.tsv.gz`.  
   - Reads the counts table with `data.table::fread`.  
   - Uses `GeneID` as row names and removes it from the count matrix.

3. **Selects samples**  
   - Removes intermediate‑passage samples for MRC‑5 and HFF, keeping young and old samples as in the original design of the study.

4. **Filters lowly expressed genes**  
   - Removes genes with zero counts in all samples.  
   - Keeps genes with at least a minimum CPM in a minimum fraction of samples (default: CPM ≥ 3 in ≥ 25% of samples).
   - Returns a filtered counts matrix `GSE63577_filtered`.

5. **Generates basic QC plots**  
   Writes plots into `filtering_result/`:
   - Density plots of log‑CPM before and after filtering.  
   - Histograms of “number of samples per gene” and “number of expressed genes per sample”.  
   - A log–log scatter plot of library size vs number of expressed genes.

6. **Creates sample and gene annotations**  
   - `sample_info`: `SampleID`, `Cell_type`, `Passage` (young/old) for each GSM sample.  
   - `gene_info`: `GeneID` for all filtered genes.  
   - Row names are set appropriately for downstream use with Bioconductor classes such as `SummarizedExperiment`.

**Main objects created**

- `GSE63577_filtered` – filtered gene × sample counts matrix.  
- `sample_info` – sample annotation data frame.  
- `gene_info` – gene annotation data frame.

---

## Script 2 – Normalization, QC and differential expression (`GSE63577_normalization_DEG.R`)

This script takes the filtered counts (`GSE63577_filtered`), `sample_info` and `gene_info` from Script 1 and performs normalization and DEG analysis.

### Normalization and QC (`normalize_counts()`)

The function `normalize_counts()`:

- Checks that `rownames(X)` match `rownames(gene.annotation)` and `colnames(X)` match `rownames(sample.annotation)`.  
- Builds:
  - An `edgeR::DGEList` for TMM normalization and log‑CPM.  
  - A `DESeq2::DESeqDataSet` for size factor and dispersion estimation and for variance‑stabilizing transformation (VST).
- Applies TMM normalization (`calcNormFactors(method = "TMM")`) on the `DGEList` and computes log‑CPM values (`cpm(log = TRUE)`).
- Runs `DESeq2::vst` (blind = TRUE and optionally blind = FALSE) on the DESeq2 object to obtain variance‑stabilized data for exploratory analysis.
- Produces QC plots:
  - mean–SD plots for log‑CPM and vst values.  
  - Library size barplot.  
  - RLE plots for log‑CPM and vst (Relative Log Expression).  
  - PCA (on log‑CPM and vst) and MDS plots to check sample clustering and potential batch effects.

**Main outputs from `normalize_counts()`**

- `dge` – `edgeR::DGEList` with TMM factors.  
- `lcpm` – TMM‑normalized log‑CPM matrix.  
- `dds` – `DESeq2::DESeqDataSet`.  
- `vsd`, `vsdbf` – vst‑transformed data (blind TRUE/FALSE).  
- `res_pca`, `res_pca_vst` – PCA results on log‑CPM and vst values.  
- `var_feat`, `var_feat_vst` – selected features for PCA.

QC plots are saved under `normalizing_result/`.

### Differential expression

DE analysis uses `omicsTools::differential_expression()` built on edgeR/limma.

- A design matrix is built as `~ 0 + Passage`, so each cell type × passage combination (e.g. BJyoung, BJold, etc.) gets a dedicated column.
- A contrast matrix is created with `limma::makeContrasts` to compare old vs young for each fibroblast line:

  - HFF: old vs young  
  - MRC‑5: old vs young  
  - IMR‑90: old vs young  
  - BJ: old vs young  
  - WI‑38: old vs young

- `omicsTools::differential_expression()` is run with:
  - `dge` – TMM‑normalized `DGEList`  
  - `design` – design matrix  
  - `contr_mat` – contrast matrix  

Results (fit object and per‑contrast top tables) are saved under `DEG_result/` and in the R object `df_deg`.

**Main objects created**

- `norm_res` – list returned by `normalize_counts()`.  
- `dge` – `DGEList` used for DE.  
- `df_deg` – list containing:
  - `fit` – limma fit object.  
  - `tt` – per‑contrast top tables (e.g. `df_deg$tt$oHFF_yHFF`).

---

## Script 3 – Post‑processing and visualization (`GSE63577_postprocessing_visualization.R`)

This script operates on DE results and raw counts to create annotated tables and visual summaries.

### Annotated raw matrix

- Downloads the same raw count matrix (`GSE63577_raw_counts_GRCh38.p13_NCBI.tsv.gz`) as in Script 1 and removes intermediate samples.  
- Merges it with an external annotation table `Homo_sapiens.GRCh38.110.gtf.entrez.xlsx` containing fields such as `gene_id`, `gene_name`, `gene_biotype`, and Entrez IDs.
- Exports the fully annotated raw counts as `results_to_do_list/Raw_Matrix_GSE63577_annotated.xlsx`.

### DEG tables with flags

- For each contrast (HFF, MRC‑5, IMR‑90, BJ, WI‑38) it:
  - Extracts DE results from `df_deg`.  
  - Merges them with gene annotation.  
  - Classifies genes as `up`, `down`, or `NS` using thresholds `adj.P.Val < 0.05` and \(|\text{logFC}| > 0.58\).
  - Exports the annotated DEG table (including flags) as Excel.

### lncRNA and protein‑coding subsets and UpSet plots

- From each DEG table, defines subsets:
  - lncRNA (`gene_biotype == "lncRNA"`).  
  - Protein‑coding genes (`gene_biotype == "protein_coding"`). 
- For lncRNA and for protein‑coding separately:
  - Uses only `up`/`down` genes to construct gene sets per cell line.  
  - Creates a binary matrix (gene × cell line) and uses `UpSetR` to visualize overlaps of modulated genes across the five fibroblast lines (HFF, MRC‑5, IMR‑90, BJ, WI‑38).
  - Identifies:
    - Genes shared by all 5 comparisons.  
    - Genes shared by at least 3 comparisons.  
  - Merges their per‑contrast information into fully annotated tables and exports them as Excel files.

### Volcano plots

- For each contrast, uses `EnhancedVolcano` to visualize DEGs:
  - x‑axis: logFC.  
  - y‑axis: \(-\log_{10}(\text{adj.P.Val})\).  
  - Highlights genes passing `adj.P.Val < 0.05` and \(|\text{logFC}| > 0.58\).
- The title of each volcano plot reports the number of up‑regulated, down‑regulated and total significant genes.

### Heatmaps

Using the raw counts (reloaded and matched to annotation):

- **lncRNA heatmaps**:  
  - For each cell line, takes modulated lncRNAs (up/down), extracts their raw counts, log‑transforms and z‑scores across genes.  
  - Builds heatmaps (rows = lncRNAs, columns = young/old replicates) with `ComplexHeatmap::Heatmap` to show transcriptional changes in non‑coding RNAs.

- **Protein‑coding heatmaps**:  
  - For each cell line, selects the top 50 up‑regulated and top 50 down‑regulated protein‑coding genes (by logFC).
  - Extracts expression from the annotated raw counts, z‑scores per gene and generates heatmaps to highlight key senescence‑associated signatures.

**Main outputs**

- Annotated Excel tables (raw matrix, per‑contrast DEGs, lncRNA/protein‑coding overlaps).  
- UpSet plots, volcano plots, and heatmaps summarizing transcriptional changes.

---

## Script 4 – Functional enrichment (`GSE63577_enrichment.R`)

This script performs GO and KEGG enrichment for each cell line using **g:Profiler** via `gprofiler2::gost`, and summarizes pathways across fibroblast lines.

### Per‑cell line enrichment

For each contrast (`oHFF_yHFF`, `oMRC5_yMRC5`, `oIMR90_yIMR90`, `oBJ_yBJ`, `oWI38_yWI38`):

1. **Rebuilds DEG tables with flags**  
   - Uses `limma::topTable` on `fit`.  
   - Merges with annotation (`Homo_sapiens.GRCh38.110.gtf.entrez.xlsx`).  
   - Flags genes as `up`, `down`, or `NS` using the same thresholds (`adj.P.Val < 0.05`, \(|\text{logFC}| > 0.58\)).
2. **Combined enrichment (up + down)**  
   - Creates a combined list of significant DEGs (up + down), typically using Ensembl IDs (`ens_gene_id`) as the query for `gprofiler2::gost`.
   - Runs `gost` with:
     - `organism = "hsapiens"`  
     - `sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG")`  
     - `evcodes = TRUE`  
   - Saves full results as `Enrich_tbl_<CELL>_GSE63577.xlsx` in `Enrich_gprofiler/No_Dis/`.  
   - Filters KEGG terms (source == "KEGG", term_size ≤ 300), sorts by intersection size and exports `KEGG_<CELL>_GSE63577.xlsx`.  
   - Produces KEGG barplots (intersection size vs KEGG term, colored by p‑value).

3. **Separated enrichment for up‑ vs down‑regulated genes**

   - Splits DEGs into `up` and `down` based on logFC and adjusted p‑value, and sorts within each group.  
   - Uses gene symbols (`gene_name`) as queries and runs `gost` with:

     ```r
     gost(list("up-regulated" = up_gene_names,
               "down-regulated" = down_gene_names),
          organism = "hsapiens",
          sources  = c("GO:BP","GO:MF","GO:CC","KEGG"),
          evcodes  = TRUE)
     ```

   - Saves the full result table (`Enrich_tbl_o<CELL>_vs_y<CELL>_GSE63577.xlsx`) under `Enrich_gprofiler/`.  
   - Extracts KEGG results into `KEGG_o<CELL>_y<CELL>_GSE63577.xlsx`.  
   - Creates faceted KEGG barplots (separate panels for up‑ and down‑regulated genes), with intersection size on the x‑axis and KEGG terms on the y‑axis.

### Combined KEGG dotplots across fibroblast lines

- Collects KEGG enrichment results for all 5 cell lines:

  - From multi‑query (up vs down) results:  
    - Combines `kegg_df_<CELL>` tables into one data frame with a `cell_line` and `regulation` column.  
    - Writes `KEGG_5Cell_type_up_down_GSE63577.xlsx`.  
    - Plots a dotplot where:
      - x‑axis: cell line (BJ, IMR90, WI38, MRC5, HFF).  
      - y‑axis: KEGG term.  
      - Dot size: intersection size.  
      - Dot color: regulation (up‑regulated vs down‑regulated).

  - From single‑query combined results (no up/down separation):  
    - Combines KEGG tables from combined analyses into `KEGG_No_Dis_5Cell_type_GSE63577.xlsx` in `Enrich_gprofiler/No_Dis/`.  
    - Draws a dotplot where dot color reflects p‑value (gradient) instead of regulation, useful when only overall enrichment is of interest.

**Main outputs**

- Per‑cell line enrichment tables for combined and separated (up/down) gene sets (GO and KEGG).  
- KEGG barplots per cell line.  
- Combined KEGG dotplots summarizing pathway usage across fibroblast lines.

---

## Requirements

You will need a recent version of R and the following packages:

**CRAN**

- `data.table`, `readxl`, `openxlsx`, `dplyr`, `ggplot2`

**Bioconductor / others**

- `edgeR`, `DESeq2`, `SummarizedExperiment`, `Biobase`, `BiocGenerics`  
- `pcaMethods`, `plotrix`, `vsn`, `pals`, `ComplexHeatmap`, `limma`  
- `omicsTools` (for DE helper functions)  
- `UpSetR`, `EnhancedVolcano`, `gprofiler2`

Example installation:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "edgeR",
  "DESeq2",
  "SummarizedExperiment",
  "Biobase",
  "BiocGenerics",
  "pcaMethods",
  "plotrix",
  "vsn",
  "pals",
  "ComplexHeatmap",
  "limma",
  "omicsTools"
))

install.packages(c(
  "data.table",
  "readxl",
  "openxlsx",
  "dplyr",
  "ggplot2",
  "UpSetR",
  "EnhancedVolcano",
  "gprofiler2"
))
```

---

## How to run

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/<your-repo-name>.git
cd <your-repo-name>
```

### 2. Preprocessing (Script 1)

```r
source("GSE63577_preprocessing.R")

# Created:
GSE63577_filtered  # filtered counts
sample_info        # sample annotation
gene_info          # gene annotation
```

### 3. Normalization and DE (Script 2)

```r
source("GSE63577_normalization_DEG.R")

# Created:
norm_res  # list from normalize_counts()
dge       # edgeR DGEList for DE
df_deg    # DE results (fit + per-contrast tables)
```

### 4. Post‑processing and visualization (Script 3)

```r
source("GSE63577_postprocessing_visualization.R")
```

Outputs include:

- Annotated raw matrix (`Raw_Matrix_GSE63577_annotated.xlsx`).  
- Per‑contrast DEG tables with flags.  
- lncRNA/protein‑coding overlap tables (shared across cell lines).  
- UpSet plots, volcano plots and heatmaps for modulated lncRNAs and protein‑coding genes.

### 5. Functional enrichment (Script 4)

```r
source("GSE63577_enrichment.R")
```

Outputs include:

- GO/KEGG enrichment tables per cell line (combined and up/down separated).  
- KEGG barplots per fibroblast line.  
- Combined KEGG dotplots summarizing enriched pathways across all five cell lines.

---

## Suggested repository structure

```text
.
├── GSE63577_preprocessing.R
├── GSE63577_normalization_DEG.R
├── GSE63577_postprocessing_visualization.R
├── GSE63577_enrichment.R
├── Homo_sapiens.GRCh38.110.gtf.entrez.xlsx
├── README.md
├── filtering_result/        # QC plots from Script 1
├── normalizing_result/      # QC plots from Script 2
├── DEG_result/              # DEG outputs from Script 2
├── results_to_do_list/      # post-processing outputs from Script 3
└── Enrich_gprofiler/        # enrichment outputs from Script 4 (and Enrich_gprofiler/No_Dis/)
```

---

## Citation

If you use this repository or scripts in a publication, please cite:

- The GEO dataset: **GSE63577 – RNA‑seq of human fibroblasts during replicative senescence**.
- The original study describing conserved senescence‑associated genes and pathways in primary human fibroblasts.
- g:Profiler for enrichment analysis (`gprofiler2` / g:Profiler toolset).
