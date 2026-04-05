# 🧬 GSE63577 RNA-seq Analysis (Human Fibroblast Senescence)

This repository contains a complete RNA-seq analysis pipeline for the dataset **GSE63577**, focusing on transcriptional changes during replicative senescence in human fibroblasts.

The workflow includes preprocessing, normalization, differential expression analysis, visualization, and functional enrichment.

---

## 📊 Dataset

- **Accession:** GSE63577  
- **Source:** GEO (NCBI)  
- **Organism:** *Homo sapiens*  
- **Cell lines:** BJ, WI-38, IMR-90, HFF, MRC-5  
- **Conditions:** Young vs Senescent (replicative aging)

---

## ⚙️ Pipeline Overview

The analysis is organized into four main scripts:

### 1️⃣ Preprocessing
**`GSE63577_preprocessing.R`**

- Download raw counts from GEO  
- Remove unwanted samples  
- Filter lowly expressed genes (CPM-based)  
- Generate QC plots  
- Create:
  - `GSE63577_filtered`
  - `sample_info`
  - `gene_info`

---

### 2️⃣ Normalization & Differential Expression
**`GSE63577_normalization_DEG.R`**

- TMM normalization (edgeR)  
- VST transformation (DESeq2)  
- QC:
  - PCA
  - MDS
  - RLE
- Differential expression using limma/edgeR  
- Outputs:
  - `dge`
  - `df_deg` (DE results)

---

### 3️⃣ Post-processing & Visualization
**`GSE63577_postprocessing_visualization.R`**

- Annotate genes (gene name, biotype, Entrez)  
- Classify DEGs:
  - up / down / NS  
- Generate:
  - Volcano plots  
  - Heatmaps  
  - UpSet plots  
- Separate analysis for:
  - lncRNAs  
  - Protein-coding genes  

---

### 4️⃣ Functional Enrichment
**`GSE63577_enrichment.R`**

- GO & KEGG enrichment using `gprofiler2`  
- Analysis:
  - Combined DEGs  
  - Up vs Down separately  
- Outputs:
  - Enrichment tables (Excel)  
  - KEGG barplots  
  - Cross-cell-line dotplots  

---

## 📁 Output Structure
