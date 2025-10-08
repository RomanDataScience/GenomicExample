#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import statsmodels.api as sm

INPUT = "data/raw/mock_allen.h5ad"
OUTDIR = Path("results")
OUTDIR.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(INPUT)
print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[adata.obs["n_genes_by_counts"] > 200, :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1000)
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)

sc.pl.umap(adata, color=["cell_type", "diagnosis"], save="_mock_pipeline.png", show=False)

adata.obs["group"] = adata.obs[["donor_id","cell_type","diagnosis"]].astype(str).agg("_".join, axis=1)
pb_counts = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
pb_grouped = pb_counts.groupby(adata.obs["group"]).sum().T
meta = adata.obs.drop_duplicates("group").set_index("group")[["cell_type","diagnosis","age","sex"]]

pb_grouped.to_csv("results/tables/pseudobulk_counts.csv")
meta.to_csv("results/tables/pseudobulk_metadata.csv")

def run_de(gene_expr, meta):
    meta["diagnosis_bin"] = (meta["diagnosis"] == "AD").astype(int)
    meta["sex_bin"] = (meta["sex"] == "M").astype(int)
    X = sm.add_constant(meta[["diagnosis_bin", "age", "sex_bin"]])
    model = sm.OLS(gene_expr, X).fit()
    return model.params["diagnosis_bin"], model.pvalues["diagnosis_bin"]

results = []
for gene in pb_grouped.index[:200]:
    beta, pval = run_de(pb_grouped.loc[gene], meta)
    results.append({"gene": gene, "beta_diagnosis": beta, "pval": pval})
de = pd.DataFrame(results).sort_values("pval")
de.to_csv("results/tables/de_results_mock.csv", index=False)

print("✅ Analysis complete. Results saved in results/tables/")
