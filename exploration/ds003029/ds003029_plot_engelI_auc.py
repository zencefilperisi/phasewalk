"""
ds003029 - figure: per-subject CTQW AUC in the Engel-I group (the primary
result of ADR 0013). Reads cohort_results.csv (produced by the cohort analysis)
and saves phase3_engelI_auc.png. Points below the chance line (0.5) are red.

Requires: pandas, numpy, matplotlib
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
df = pd.read_csv(HERE / "cohort_results.csv")

e1 = (df[df["engel"] == 1]
      .dropna(subset=["auc_ctqw_pr"])
      .sort_values("auc_ctqw_pr")
      .reset_index(drop=True))
auc = e1["auc_ctqw_pr"].values
med = float(np.median(auc))

fig, ax = plt.subplots(figsize=(8, 5.2))
x = np.arange(len(auc))
colors = ["#C44E52" if a < 0.5 else "#4C72B0" for a in auc]
ax.scatter(x, auc, s=95, c=colors, zorder=3, edgecolor="white", linewidth=0.9)
ax.axhline(0.5, color="#555", ls="--", lw=1.3, zorder=1)
ax.text(len(auc) - 0.5, 0.515, "chance (0.5)", va="bottom", ha="right",
        fontsize=9, color="#555")
ax.axhline(med, color="#C44E52", ls=":", lw=1.4, zorder=1)
ax.text(0.2, med + 0.03, f"median {med:.2f}", va="bottom", ha="left",
        fontsize=9, color="#C44E52")
ax.set_xticks(x)
ax.set_xticklabels(e1["label"], rotation=90, fontsize=7.5)
ax.set_ylim(0, 1)
ax.set_ylabel("AUC  (low PR → SOZ)", fontsize=10)
ax.set_title("CTQW does not localize the seizure onset zone above chance\n"
             "Engel I (seizure-free) patients, n=17", fontsize=11.5, pad=12)
ax.spines[["top", "right"]].set_visible(False)
ax.margins(x=0.04)
plt.tight_layout()
out = HERE / "phase3_engelI_auc.png"
plt.savefig(out, dpi=170, bbox_inches="tight")
print(f"saved -> {out}  (median {med:.3f}, {(auc < 0.5).sum()}/{len(auc)} below chance)")