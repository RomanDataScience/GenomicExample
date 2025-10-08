# 🧠 Allen Brain Alzheimer’s Transcriptomic Analysis  
**Cell-Type–Specific Differential Expression and Pathway Discovery Using Allen Brain (SEA-AD) Data**

---

## 📘 Overview

This repository provides a **reproducible research pipeline** for analyzing **single-cell and single-nucleus transcriptomic data** from the **Allen Brain Institute’s SEA-AD dataset** (or synthetic mock data) to identify **cell-type–specific molecular signatures of Alzheimer’s disease**.

The project demonstrates a full bioinformatics workflow — from **data ingestion and QC** to **dimensionality reduction**, **pseudo-bulk differential expression (DE)**, and **functional enrichment analysis** — using open-source tools in Python and R.

It includes:
- 🧬 A **synthetic dataset generator** to simulate Allen-like data.
- 🧹 Scripts for **QC, normalization, and UMAP visualization**.
- 📊 **Pseudo-bulk aggregation** by donor and cell type.
- 🧮 **Covariate-adjusted DE** analysis (diagnosis + age + sex).
- 📈 **Functional enrichment** via g:Profiler.
- 🐳 A **Dockerized environment** for full reproducibility.

---

## 🏗️ Project Structure

allen-brain-ad-analysis/
├── README.md
├── environment/
│ ├── Dockerfile
│ ├── requirements.txt
│ └── env.yml
├── data/
│ ├── raw/ # Raw or synthetic data (mock_allen.h5ad)
│ └── processed/ # QC’d and processed datasets
├── notebooks/
│ ├── 00_generate_mock_allen_data.ipynb
│ ├── 01_qc_and_preprocessing.ipynb
│ ├── 02_umap_and_clustering.ipynb
│ ├── 03_pseudobulk_DE.ipynb
│ └── 04_functional_analysis.ipynb
├── scripts/
│ ├── allen_pipeline.py # Main end-to-end workflow
│ ├── run_pseudobulk_DE.R # R DESeq2 analysis
│ ├── pathway_enrichment.py # Pathway enrichment
│ ├── utils.py # Helper functions
│ └── config.yaml # Config parameters
├── results/
│ ├── figures/
│ ├── tables/
│ └── logs/
└── .gitignore