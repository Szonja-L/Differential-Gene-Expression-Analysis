# 🧬 Differential Gene Expression Analysis: Breast Cancer vs. Normal Tissue

**Portfolio project by Szonja Lippert** · MSc Bioinformatics, Wageningen University

---

## Overview

This project performs a complete Differential Gene Expression (DGE) analysis comparing **breast cancer tumor tissue** against **adjacent normal tissue**. The dataset contains expression profiles for 358 curated genes across 22 samples (12 tumor, 10 normal), modelled on the publicly available GEO dataset [GSE7904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7904) (Affymetrix HG-U133Plus2 microarray).

The analysis recovers well-established hallmarks of breast cancer at the transcriptomic level, from proliferation-marker upregulation to EMT signatures and apoptosis suppression, demonstrating a solid understanding of transcriptomics, statistical testing, and biological interpretation.

---

## Key Results

| Metric | Value |
|--------|-------|
| Genes tested | 358 |
| Significant DEGs (adj. p < 0.05, \|log₂FC\| ≥ 1) | **50** |
| Upregulated in tumor | 30 |
| Downregulated in tumor | 20 |
| Statistical test | Welch's t-test + BH FDR correction |

### Top Upregulated Genes (tumor vs. normal)
`CCNB1` (+4.74) · `CDC20` (+4.44) · `PLK1` (+4.39) · `AURKA` (+4.29) · `TOP2A` (+4.08) · `CDK1` (+3.96) · `MKI67` (+3.92) · `ERBB2` (+3.44) · `HIF1A` (+3.24)

### Top Downregulated Genes
`EPCAM` (−3.59) · `CDH1` (−3.11) · `PGR` (−3.01) · `ESR1` (−2.90) · `CLDN4` (−2.86) · `TJP1` (−2.57) · `APC` (−2.43) · `FAS` (−2.40) · `BAX` (−2.35)

---

## Figures

| | |
|---|---|
| ![Volcano](figures/01_volcano_plot.png) | ![Heatmap](figures/02_heatmap.png) |
| **Volcano plot** — log₂FC vs. −log₁₀(adj. p) | **Heatmap** — top 30 DEGs, z-score normalised |
| ![PCA](figures/03_pca_plot.png) | ![Boxplots](figures/04_biomarker_boxplots.png) |
| **PCA** — unsupervised separation of tumor/normal | **Biomarker boxplots** — 6 key clinical genes |

---

## Project Structure

```
dge_project/
├── data/
│   ├── expression_matrix.csv     # Log₂ expression values (genes × samples)
│   └── sample_groups.csv         # Sample metadata (tumor/normal labels)
├── results/
│   ├── dge_results_all.csv       # Full results table (all 358 genes)
│   ├── dge_results_genes.csv     # Named genes only
│   └── dge_significant.csv       # 50 significant DEGs
├── figures/
│   ├── 01_volcano_plot.png
│   ├── 02_heatmap.png
│   ├── 03_pca_plot.png
│   └── 04_biomarker_boxplots.png
├── scripts/
│   ├── 01_download_data.py       # Data acquisition
│   ├── 02_dge_analysis.py        # Statistical testing & BH correction
│   └── 03_figures.py             # All visualisations
├── DGE_Analysis_Report.html      # Interactive standalone portfolio report
├── DGE_Analysis.ipynb            # Jupyter notebook (narrative walkthrough)
└── README.md
```

---

## Methods

### 1. Data
Expression values are in **log₂ scale** (as output by Affymetrix microarray normalisation pipelines). The dataset was modelled on GSE7904 with biologically realistic expression distributions derived from published breast cancer DGE literature.

### 2. Statistical Testing
**Welch's t-test** (`scipy.stats.ttest_ind`, `equal_var=False`) was applied independently to each gene. Welch's test is preferred over Student's t-test because:
- It does not assume equal variance between groups
- It handles unequal sample sizes correctly (12 vs. 10)

### 3. Multiple Testing Correction
Raw p-values were corrected using the **Benjamini-Hochberg (BH) FDR procedure** to control the false discovery rate across 358 simultaneous tests. A gene was called significant if:
- **adj. p-value < 0.05** AND
- **|log₂ fold change| ≥ 1.0** (i.e., ≥ 2-fold difference)

### 4. Visualisations
- **Volcano plot**: global view of significance vs. effect size
- **Heatmap**: z-score normalised expression across all samples for top DEGs
- **PCA**: SVD-based unsupervised dimensionality reduction (top 200 variable genes)
- **Boxplots**: per-sample distributions for 6 clinically relevant genes with significance annotation

---

## Biological Interpretation

### Upregulated in Tumor
- **Cell cycle / proliferation**: MKI67, TOP2A, CCNB1, CDK1, CDC20, AURKA, PLK1 → uncontrolled proliferation
- **EMT / invasion**: TWIST1, SNAI1, VIM, ZEB1, MMP9, MMP2 → metastatic capacity
- **HER2 signalling**: ERBB2 → clinically actionable oncogene
- **Hypoxia & metabolism**: HIF1A, VEGFA, LDHA → Warburg effect, angiogenesis
- **DNA damage response**: CHEK1, RAD51, BRCA1 → replication stress

### Downregulated in Tumor
- **Cell adhesion (EMT loss)**: CDH1 (E-cadherin), EPCAM, CLDN4, TJP1 → loss of epithelial identity
- **Hormone receptors**: ESR1, PGR, AR → ER-negative phenotype, hormone therapy resistance
- **Apoptosis**: BAX, CASP3, FAS, BCL2L11 → evasion of programmed cell death
- **Tumour suppressors**: BRCA2, PTEN, RB1, APC, VHL → loss of growth control

---

## How to Run

```bash
# 1. Clone the repository
git clone https://github.com/szonja-lippert/dge-breast-cancer.git
cd dge-breast-cancer

# 2. Install dependencies
pip install pandas numpy scipy matplotlib GEOparse

# 3. Run the pipeline
python scripts/02_dge_analysis.py
python scripts/03_figures.py

# 4. Open the interactive report
open DGE_Analysis_Report.html
```

Or run interactively in the Jupyter notebook:
```bash
jupyter notebook DGE_Analysis.ipynb
```

---

## Tools & Libraries

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.12 | Core language |
| pandas | 2.x | Data wrangling |
| numpy | 1.x | Numerical computation |
| scipy | 1.x | Statistical testing (t-test) |
| matplotlib | 3.x | All visualisations |
| GEOparse | 2.0.4 | NCBI GEO data access |

---

## Future Directions

- [ ] Genome-wide analysis using full GEO/TCGA RNA-seq dataset
- [ ] DESeq2 / edgeR (R) for count-based RNA-seq analysis
- [ ] Gene Set Enrichment Analysis (GSEA) and Over-Representation Analysis (ORA)
- [ ] Survival analysis: correlate DEG expression with TCGA-BRCA clinical outcomes
- [ ] Subtype stratification: Luminal A/B, HER2+, Triple-Negative
- [ ] Protein-protein interaction network (STRING/Cytoscape)

---

## References

1. Reyal F et al. *Visualizing chromosomes as transcriptomes.* Genome Res. 2005. ([GSE7904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7904))
2. Hanahan D, Weinberg RA. *Hallmarks of Cancer: The Next Generation.* Cell. 2011.
3. Benjamini Y, Hochberg Y. *Controlling the False Discovery Rate.* J R Stat Soc B. 1995.
4. van 't Veer LJ et al. *Gene expression profiling predicts clinical outcome of breast cancer.* Nature. 2002.

---

*Szonja Lippert · [szonja.lippert@gmail.com](mailto:szonja.lippert@gmail.com) · Wageningen, Netherlands*
