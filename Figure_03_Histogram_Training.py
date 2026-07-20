#!/usr/bin/env python3
"""
Figure 3 (left) — Binding site length distribution in training datasets

Categories:
- Very short (≤30 AA)
- Medium short (31–70 AA)
- Long (>70 AA)

Author: Nawar Malhis
Refined: Grok
The University of British Columbia, 2026
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# Add AFF project path
from param import *
if aff_path not in globals() or aff_path not in sys.path:  # type: ignore[name-defined]
    sys.path.append(str(aff_path))

from annotated_fasta import aff_load3, aff_tag_size, aff_gen_counts


def load_binding_site_fractions() -> dict:
    """Compute percentage of binding residues in each length category."""
    data_config = {
        "DisProt": {"tag": "binding_protein", "size_ranges": {
            "Very_short": [0, 31],
            "Medium_short": [31, 71],
            "Long": [71, 100000]
        }},
        "CAID1u": {"tag": "binding_protein", "size_ranges": {
            "Very_short": [0, 31],
            "Medium_short": [31, 71],
            "Long": [71, 100000]
        }},
        "CAID23u": {"tag": "binding_protein", "size_ranges": {
            "Very_short": [0, 31],
            "Medium_short": [31, 71],
            "Long": [71, 100000]
        }},
        "DBs": {"tag": "PDB", "size_ranges": {
            "Very_short": [0, 31],
            "Medium_short": [31, 71],
            "Long": [71, 100000]
        }},
    }

    fractions = {"Very_short": [], "Medium_short": [], "Long": []}
    datasets = ["DisProt", "CAID1u", "CAID23u", "DBs", "TR2008u", "TR2017"]

    for in_data in ["DisProt", "CAID1u", "CAID23u", "DBs"]:
        cfg = data_config[in_data]
        tag = cfg["tag"]

        # Load once to get total binding residues
        af = aff_load3(f"Data/af/{in_data}.af")
        aff_gen_counts(af)
        total = af["metadata"]["counts"]["tags_dict"][tag]["1"]

        for category, sz_range in cfg["size_ranges"].items():
            af = aff_load3(f"Data/af/{in_data}.af")
            aff_tag_size(af, tag, sz_range=sz_range)
            aff_gen_counts(af)
            cnt = af["metadata"]["counts"]["tags_dict"][tag]["1"]
            percentage = round(cnt * 100.0 / total, 1) if total > 0 else 0.0
            fractions[category].append(percentage)

    # TR2008u & TR2017: 100% very short (as per preprint)
    fractions["Very_short"].extend([100.0, 100.0])
    fractions["Medium_short"].extend([0.0, 0.0])
    fractions["Long"].extend([0.0, 0.0])

    return fractions, datasets


if __name__ == "__main__":
    fractions, datasets = load_binding_site_fractions()

    # Validation
    print("Binding residue length fractions (%):")
    for i, ds in enumerate(datasets):
        total = sum(fractions[cat][i] for cat in fractions)
        print(f"{ds:10} → Very_short: {fractions['Very_short'][i]:5.1f} | "
              f"Medium: {fractions['Medium_short'][i]:5.1f} | "
              f"Long: {fractions['Long'][i]:5.1f}  (sum={total:.1f}%)")

    # ========================== PLOTTING ==========================
    fig, ax = plt.subplots(figsize=(10, 6.5))

    x = np.arange(len(datasets))
    width = 0.78

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

    # Stacked bars
    bottom = np.zeros(len(datasets))
    for i, (cat, color) in enumerate(zip(fractions.keys(), colors)):
        values = fractions[cat]
        ax.bar(x, values, width, bottom=bottom, label=cat.replace("_", " "),
               color=color, edgecolor="white", linewidth=0.5)
        bottom += values

    # Formatting
    ax.set_ylabel("% of binding residues", fontweight="bold", fontsize=15)
    ax.set_title("Binding site length distribution: Training datasets",
                 fontsize=16, pad=25)

    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=14, rotation=0)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_ylim(0, 105)

    ax.legend(loc="upper center", fontsize=13, ncol=3, bbox_to_anchor=(0.5, -0.08))
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    # Percentage labels inside bars
    for i in range(len(datasets)):
        cum = 0.0
        for cat in fractions:
            perc = fractions[cat][i]
            if perc > 4.0:   # label only larger segments
                ax.text(i, cum + perc / 2, f"{perc}%", ha="center", va="center",
                        fontsize=11.5, color="white", fontweight="bold")
            cum += perc

    plt.tight_layout()

    # Save
    out_path = "Data/results/Figure_3/Figure_3_left_training.png"
    plt.savefig(out_path, dpi=400, bbox_inches="tight")
    print(f"Figure saved: {out_path}")

    plt.show()
