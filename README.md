# README — Genomic Data Analysis Case Study Using Allen Brain Data

## Project: Discovering Cell-Type Specific Gene Expression Changes in Alzheimer’s Disease

### Data Source
Allen Brain Atlas / Allen Brain Institute Cell Types Database  
Transcriptomic and anatomical/geospatial gene expression data  
Resources: [Allen Brain Map](https://portal.brain-map.org/), PMC references, Cell Types Database.

---

### Goals
- Identify genes whose expression patterns are altered in specific brain regions or cell types in Alzheimer’s disease (AD) versus controls.  
- Understand which biological pathways are implicated, particularly cell-type specific contributions (neuronal, glial, etc.).  
- Provide hypotheses for functional follow-ups (e.g., potential drug targets or biomarkers).

---

### Steps / Methods

#### 1. Data Preparation & QC
- Obtain transcriptomic datasets from the Allen Brain Institute for human brain samples (healthy and AD, or matched anatomical references).
- Preprocess data: normalization, filtering of low-expression genes, correction for batch effects, and alignment of metadata (age, sex, brain region).

#### 2. Exploratory Data Analysis
- Dimensionality reduction (PCA / UMAP) to visualize clustering by brain region, disease status, or cell type.  
- Hierarchical clustering of regions/cell types using expression profiles.  
- Identify differentially expressed genes (DEGs) between AD vs. controls in selected brain regions or cell types.

#### 3. Cell-Type Deconvolution or Spatial Mapping
- Use cell-type marker genes (from Allen Cell Types data) to infer cell-type proportions in bulk RNA-seq samples.  
- Alternatively, leverage single-cell or single-nucleus RNA-seq datasets from Allen for cell-type specific analysis.  
- Map differential expression patterns to specific cell types: e.g., upregulation of microglial genes or downregulation of excitatory neuronal synaptic genes.

#### 4. Functional Interpretation
- Perform enrichment analysis on DEGs (Gene Ontology, KEGG pathways, Reactome).  
- Investigate spatial patterns of altered genes—do they cluster in regions like hippocampus or cortex?  
- Construct co-expression networks to identify modules or hub genes disrupted in AD.

#### 5. Custom Modeling / Visualization
- Build interactive heatmaps and spatial expression maps overlayed on anatomical brain regions.  
- Explore machine learning models to predict disease status from integrated cell-type expression patterns.  
- Develop visualization tools to allow collaborators to explore results by region or cell type.

#### 6. Collaborative Guidance / Transferability
- Train team members to use and extend these analysis pipelines.  
- Provide reproducible workflows (Jupyter notebooks, version-controlled scripts, and Docker/conda environments).  
- Document all steps for transparency and reproducibility.

---

### Expected Results / Insights
- Discovery of genes differentially expressed in AD, enriched in processes like microglial activation or synaptic dysfunction.  
- Identification of key cell types driving disease progression (e.g., astrocytes, microglia).  
- Spatial maps illustrating progressive gene expression changes across cortical and hippocampal regions.  
- Network-based identification of hub genes and signaling pathways linked to neuroinflammation and neurodegeneration.  
- Creation of interactive tools and visualizations to explore expression data by region or cell type.

---

### Repository Contents
- `allen_seaad_pipeline.py`: example notebook performing preprocessing → UMAP → pseudo-bulk DE.  
- `Dockerfile` and `requirements.txt`: for full reproducibility.  
- `pseudobulk_counts.csv`, `pseudobulk_metadata.csv`, `de_results_*.csv`: example outputs.  
- `README.md` (this file): project overview and methods.

---

### Installation & Usage

#### **Option 1 — Run via Docker (Recommended)**

1. **Build the container:**
   ```bash
   docker build -t allen-seaad .
   ```

2. **Run the analysis:**
   ```bash
   docker run -v $(pwd):/app allen-seaad
   ```

3. **Outputs:**
   - Pseudo-bulk expression matrices: `pseudobulk_counts.csv`, `pseudobulk_metadata.csv`  
   - Differential expression results: `de_results_<cell_type>_covariates.csv`

---

#### **Option 2 — Run locally with Python**

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Run pipeline:**
   ```bash
   python allen_seaad_pipeline.py
   ```

3. **Visualize UMAP plots and review DE results** generated in the working directory.

---

### Citation
If you use this pipeline, please cite:
- Allen Institute for Brain Science — [SEA-AD MTG dataset](https://portal.brain-map.org/explore/transcriptome/mcg-sea-ad)
