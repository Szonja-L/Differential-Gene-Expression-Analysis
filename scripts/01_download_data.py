"""
Step 1: Download and Preprocess REAL Differential Gene Expression Data
Dataset: GSE7904 (Richardson et al.) — Affymetrix HG-U133 Plus 2 (GPL570)
Source: NCBI Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7904)

This replaces the earlier synthetic-data generator. It downloads the real
series, real sample labels, and real platform annotation directly from GEO,
then preprocesses it into the same expression_matrix.csv / sample_groups.csv
format that 02_dge_analysis.py and 03_figures.py already expect — so the
rest of the pipeline runs completely unchanged on real values.

REQUIRES INTERNET ACCESS TO ftp.ncbi.nlm.nih.gov.
Run this on your own machine (not in a sandboxed/restricted environment).
"""

import os
import numpy as np
import pandas as pd
import GEOparse

DATA_DIR = "data"
CACHE_DIR = "geo_cache"
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)

print("=" * 60)
print("DOWNLOADING GSE7904 FROM NCBI GEO")
print("=" * 60)

gse = GEOparse.get_GEO(geo="GSE7904", destdir=CACHE_DIR)

# ──────────────────────────────────────────────────────────────
# 1. Assign Tumor / Normal labels from the REAL sample titles
# ──────────────────────────────────────────────────────────────
# Real titles look like: "Basal (T118)", "BRCA1 (T151)", "Non-BLC (T183)",
# "Normal breast (B42)", "Normal organelle (NE9)".
#
# "Normal organelle" samples are excluded from the main comparison: GEO's
# own metadata for these is ambiguous about tissue origin (likely sorted
# normal epithelial cells rather than whole breast tissue), so pooling them
# with "Normal breast" would risk comparing two different cell populations
# under one label. This is a deliberate, documented judgment call — see
# README "Data & Limitations".
meta = []
for gsm_name, gsm in gse.gsms.items():
    title = gsm.metadata['title'][0]
    if "Normal breast" in title:
        group = "Normal"
    elif "Normal organelle" in title:
        group = "Exclude"
    else:
        group = "Tumor"  # Basal, BRCA1, Non-BLC subtypes
    meta.append({"gsm": gsm_name, "title": title, "group": group})

meta_df = pd.DataFrame(meta)
n_excluded = (meta_df["group"] == "Exclude").sum()
meta_df = meta_df[meta_df["group"] != "Exclude"].reset_index(drop=True)

# Friendly sample names (Tumor_01, Tumor_02, ..., Normal_01, ...) so the
# existing heatmap/labeling code in 03_figures.py keeps working unchanged.
meta_df = meta_df.sort_values(["group", "gsm"]).reset_index(drop=True)
meta_df["sample"] = (
    meta_df.groupby("group").cumcount().add(1).astype(str).str.zfill(2)
)
meta_df["sample"] = meta_df["group"] + "_" + meta_df["sample"]

print(f"\nExcluded {n_excluded} 'Normal organelle' samples (ambiguous tissue origin)")
print(meta_df["group"].value_counts().to_string())

# ──────────────────────────────────────────────────────────────
# 2. Build the raw probe-level expression matrix
# ──────────────────────────────────────────────────────────────
print("\nBuilding probe-level expression matrix...")
expr_raw = gse.pivot_samples("VALUE")              # probes x GSM accessions
expr_raw = expr_raw[meta_df["gsm"].tolist()]
expr_raw.columns = meta_df["sample"].tolist()
print(f"  {expr_raw.shape[0]} probes x {expr_raw.shape[1]} samples")

# ──────────────────────────────────────────────────────────────
# 3. Map probes -> gene symbols using the real GPL570 annotation
# ──────────────────────────────────────────────────────────────
print("\nMapping probes to gene symbols via GPL570 platform annotation...")
platform = gse.gpls["GPL570"]
ann = platform.table.copy()

gene_col_candidates = ["Gene Symbol", "GENE_SYMBOL", "gene_assignment", "Gene symbol"]
gene_col = next((c for c in gene_col_candidates if c in ann.columns), None)
if gene_col is None:
    raise ValueError(
        f"Could not find a gene-symbol column in the GPL570 annotation table. "
        f"Available columns are: {list(ann.columns)}\n"
        f"Open this list and tell me which column actually contains gene "
        f"symbols — I'll adjust gene_col_candidates above."
    )

probe2gene = ann.set_index("ID")[gene_col].dropna().astype(str)
# Some Affymetrix annotation rows pack multiple synonyms separated by '///'
probe2gene = probe2gene.str.split("///").str[0].str.strip()
probe2gene = probe2gene[probe2gene != ""]

sample_cols = meta_df["sample"].tolist()
expr_raw = expr_raw.loc[expr_raw.index.intersection(probe2gene.index)].copy()
expr_raw["gene"] = probe2gene.loc[expr_raw.index]

# ──────────────────────────────────────────────────────────────
# 4. Collapse multiple probes per gene (keep the probe with the
#    highest average expression — standard convention for arrays
#    with probe-level redundancy)
# ──────────────────────────────────────────────────────────────
print("Collapsing multiple probes per gene...")
expr_raw["_rowmean"] = expr_raw[sample_cols].mean(axis=1)
expr_raw = expr_raw.sort_values("_rowmean", ascending=False)
expr_gene = expr_raw.drop_duplicates(subset="gene", keep="first")
expr_gene = expr_gene.set_index("gene")[sample_cols]
print(f"  {expr_gene.shape[0]} unique genes after probe collapse")

# ──────────────────────────────────────────────────────────────
# 5. Floor + log2 transform
#    These values are dChip-normalised signal intensities, NOT already
#    log2-scaled (unlike RMA/MAS5-log output). dChip's background
#    correction can also produce small negative values for low/absent
#    probes, so we floor at 1.0 before taking log2.
# ──────────────────────────────────────────────────────────────
print("Flooring at 1.0 and log2-transforming...")
expr_gene = expr_gene.clip(lower=1.0)
expr_log2 = np.log2(expr_gene)

# ──────────────────────────────────────────────────────────────
# 6. Save in the format 02_dge_analysis.py / 03_figures.py expect
# ──────────────────────────────────────────────────────────────
expr_log2.to_csv(os.path.join(DATA_DIR, "expression_matrix.csv"))
meta_df[["sample", "group", "gsm", "title"]].to_csv(
    os.path.join(DATA_DIR, "sample_groups.csv"), index=False
)

print("\n" + "=" * 60)
print(f"Saved data/expression_matrix.csv  ({expr_log2.shape[0]} genes x {expr_log2.shape[1]} samples)")
print(f"Saved data/sample_groups.csv      ({len(meta_df)} samples, real GSM accessions)")
print("=" * 60)
print("\nNext: run 02_dge_analysis.py, then 03_figures.py — unchanged.")
