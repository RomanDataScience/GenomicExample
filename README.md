# Allen Brain Alzheimer’s Transcriptomic Analysis
Comprehensive analysis pipeline to identify cell-type–specific differential expression in Alzheimer’s disease using the Allen Brain Institute’s SEA-AD dataset.

## Overview
This project reproduces a genomics workflow:
1. Preprocess and QC single-cell or single-nucleus RNA-seq data.
2. Generate embeddings (PCA, UMAP) and cluster by cell type.
3. Create pseudo-bulk matrices per donor and cell type.
4. Run differential expression (DE) controlling for covariates (age, sex).
5. Perform pathway and enrichment analysis.
6. Visualize spatial and cell-type–specific expression patterns.

## Data
- Source: [Allen Brain Atlas — SEA-AD MTG](https://portal.brain-map.org/explore/transcriptome/mcg-sea-ad)
- Download raw `.h5ad` files manually or via the Allen SDK.

## Environment
Build Docker image:
```bash
docker build -t allen-ad-analysis -f environment/Dockerfile .
```
Run container:
```bash
docker run -v $(pwd):/app allen-ad-analysis
```
Or install locally:
```bash
pip install -r environment/requirements.txt
```
## Reproducible Workflow

| Step                  | Script / Notebook                                | Output                 |
|------------------------|--------------------------------------------------|------------------------|
| QC & Preprocessing     | `01_qc_and_preprocessing.ipynb`                 | Cleaned AnnData        |
| UMAP & Clustering      | `02_umap_and_clustering.ipynb`                  | Figures, `.h5ad`       |
| Pseudo-bulk DE         | `03_pseudobulk_DE.ipynb`, `run_pseudobulk_DE.R` | CSV tables             |
| Functional Analysis    | `04_functional_analysis.ipynb`                  | Enrichment plots       |


## Citation

If you use this repo, please cite:

- Allen Institute for Brain Science — SEA-AD MTG dataset