import pandas as pd

def summarize_de_results(df, top_n=20):
    return df.sort_values("pval_diagnosis").head(top_n)
