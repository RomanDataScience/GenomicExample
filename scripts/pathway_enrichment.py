#!/usr/bin/env python3
import pandas as pd
from gprofiler import GProfiler

de = pd.read_csv("results/tables/de_results_mock.csv")
sig_genes = de.query("pval < 0.05")["gene"].tolist()

if len(sig_genes) > 5:
    gp = GProfiler(return_dataframe=True)
    enrich = gp.profile(organism='hsapiens', query=sig_genes)
    enrich.to_csv("results/tables/pathway_enrichment.csv", index=False)
    print("✅ Pathway enrichment saved.")
else:
    print("⚠️ Not enough significant genes for enrichment.")
