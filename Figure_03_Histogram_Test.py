#!/usr/bin/env python3
"""
Figure 3 (right) — Binding site length distribution in test datasets

Categories:
- Short (≤70 AA)
- Long (>70 AA)

Author: Nawar Malhis
Refined using Grok
The University of British Columbia, 2026
"""

import os
os.environ["QT_LOGGING_RULES"] = "qt.qpa.wayland.textinput=false"
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add project path
from param import *
if aff_path not in globals() or aff_path not in sys.path:
    sys.path.append(str(aff_path))

from annotated_fasta import aff_load3, aff_tag_size, aff_gen_counts


def load_test_site_fractions() -> tuple[dict, list]:
    """Compute % of binding residues in short vs long categories for test sets."""
    data_config = {
        "CAID1uh": {"tag": "binding_protein", "size_ranges": {
            "Short": [0, 71],
            "Long": [71, 100000]
        }},
        "CAID23uh": {"tag": "binding_protein", "size_ranges": {
            "Short": [0, 71],
            "Long": [71, 100000]
        }},
        "DBsh": {"tag": "PDB", "size_ranges": {
            "Short": [0, 71],
            "Long": [71, 100000]
        }},
    }

    fractions = {"Short": [], "Long": []}
    datasets = ["CAID1uh", "CAID23uh", "DBsh"]

    for in_data in datasets:
        cfg = data_config[in_data]
        tag = cfg["tag"]

        # Load once for total
        af = aff_load3(f"Data/af/{in_data}.af")
        aff_gen_counts(af)
        total = af["metadata"]["counts"]["tags_dict"][tag]["1"]

        for category, sz_range in cfg["size_ranges"].items():
            af = aff_load3(f"Data/af/{in_data}.af")  # fresh load
            aff_tag_size(af, tag, sz_range=sz_range)
            aff_gen_counts(af)
            cnt = af["metadata"]["counts"]["tags_dict"][tag]["1"]
            percentage = round(cnt * 100.0 / total, 1) if total > 0 else 0.0
            fractions[category].append(percentage)

    return fractions, datasets


if __name__ == "__main__":
    fractions, datasets = load_test_site_fractions()

    # Validation
    print("Test datasets — binding residue length fractions (%):")
    for i, ds in enumerate(datasets):
        total = fractions["Short"][i] + fractions["Long"][i]
        print(f"{ds:10} → Short: {fractions['Short'][i]:5.1f} | "
              f"Long: {fractions['Long'][i]:5.1f}  (sum={total:.1f}%)")

    # ========================== PLOTTING ==========================
    fig, ax = plt.subplots(figsize=(9, 6.2))

    x = np.arange(len(datasets))
    width = 0.78
    colors = ["#1f77b4", "#2ca02c"]

    # Stacked bars
    bottom = np.zeros(len(datasets))
    for cat, color in zip(fractions.keys(), colors):
        values = fractions[cat]
        ax.bar(
            x,
            values,
            width,
            bottom=bottom,
            label=cat,
            color=color,
            edgecolor="white",
            linewidth=0.6,
        )
        bottom += values

    # Formatting
    ax.set_ylabel("% of binding residues", fontweight="bold", fontsize=15)
    ax.set_title("Binding site length distribution: Test datasets",
                 fontsize=16, pad=20)

    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_ylim(0, 105)

    ax.legend(loc="upper center", fontsize=13.5, ncol=2, bbox_to_anchor=(0.5, -0.08))
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    # Percentage labels
    for i in range(len(datasets)):
        cum = 0.0
        for cat in fractions:
            perc = fractions[cat][i]
            if perc > 4.0:
                ax.text(
                    i,
                    cum + perc / 2,
                    f"{perc}%",
                    ha="center",
                    va="center",
                    fontsize=12,
                    color="white",
                    fontweight="bold",
                )
            cum += perc

    plt.tight_layout()

    # Save
    out_path = "Data/results/Figure_3/Figure_3_right_test.png"
    plt.savefig(out_path, dpi=400, bbox_inches="tight")
    print(f"Figure saved: {out_path}")

    plt.show()

