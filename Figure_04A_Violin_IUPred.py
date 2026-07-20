#!/usr/bin/env python3
"""
Figure 4A (left) — Horizontal violin plot of IUPred3 disorder scores
Compares Class 0 (non-binding) and Class 1 (binding) residues across datasets.

Author: Nawar Malhis
Refined using Grok
The University of British Columbia, 2026
"""

import sys
import matplotlib.pyplot as plt

# Add AFF project path
from param import *
if aff_path not in globals() or aff_path not in sys.path:  # type: ignore[name-defined]
    sys.path.append(str(aff_path))  # type: ignore[name-defined]

from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3, aff_remove_short


def get_tags_list_scores(af: dict, tags: list, score_tag: str) -> list:
    """Extract scores for specific tag/class combinations."""
    d_lst = []
    for tag_name, class_val in tags:
        sc = []
        for entry in af["data"].values():
            tags_dict = entry.get("tags", {})
            scores = entry.get("scores", {})
            if tag_name not in tags_dict or score_tag not in scores:
                continue
            tag_array = tags_dict[tag_name]
            score_array = scores[score_tag]
            for i in range(len(tag_array)):
                if tag_array[i] == class_val:
                    sc.append(score_array[i])
        d_lst.append(sc)
    return d_lst


def violin_h_plot(
    data: list,
    labels: list,
    display_means: dict | None = None,
    title: str | None = None,
    xlabel: str = "",
    f_name: str | None = None,
    t_name: str | None = None,
    fontsize: int = 18,
):
    """Horizontal violin plot with mean lines."""
    if display_means is None:
        display_means = {}

    # Write mean values to table
    if t_name:
        with open(t_name, "w", encoding="utf-8") as fout:
            for i, lbl in enumerate(labels):
                mean_val = sum(data[i]) / len(data[i]) if data[i] else 0.0
                fout.write(f"{lbl}\t{mean_val:.3f}\n")
                print(f"{lbl}\t{mean_val:.3f}")

    plt.rcParams.update({"font.size": 14})
    fig, ax = plt.subplots(figsize=(8.5, 14))
    fig.subplots_adjust(left=0.26, right=0.96, top=0.96, bottom=0.08)

    positions = list(range(1, len(data) + 1))
    ax.set_ylim(bottom=0.5, top=8.8)
    parts = ax.violinplot(
        data,
        positions=positions,
        points=200,
        showmeans=True,
        vert=False,      # horizontal
        side="high",
    )

    # Color bodies
    for i, pc in enumerate(parts["bodies"]):
        color = display_means.get(str(i), "lightblue")
        pc.set_facecolor(color)
        pc.set_edgecolor("black")
        pc.set_alpha(0.8)

    # Mean dashed lines
    for i_str, color in display_means.items():
        i = int(i_str)
        if i < len(data):
            mean_val = sum(data[i]) / len(data[i]) if data[i] else 0.0
            ax.plot([mean_val, mean_val], [0.2, len(data) + 0.8],
                    color=color, linestyle="--", linewidth=2)

    ax.yaxis.grid(True, linestyle="--", alpha=0.7)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=fontsize - 2)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    if title:
        ax.set_title(title, fontsize=fontsize + 2, pad=20)

    if f_name:
        plt.savefig(f_name, dpi=350, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    score_tag = "IUPred3"          # change to 'LIST-S2' if needed

    # Main DisProt + extra file for PDB/IDR baseline
    af = aff_load3("Data/af/DisProt_2025_06_DBs_extra.af")
    aff_remove_short(af, cut=15)
    aff_load_caid_scores(
        af, "Data/scores/", prd_list=[score_tag],
        merged=False, remove_missing_scores=True
    )

    labels_list = ["PDB", "IDR"]
    data = get_tags_list_scores(
        af, tags=[["IDR-CAID", "0"], ["IDR-CAID", "1"]], score_tag=score_tag
    )

    # Test datasets
    d_set_dict = {
        "DBsh": "PDB",
        "CAID23uh": "binding_protein",
        "CAID1uh": "binding_protein",
    }

    for ds, tag in d_set_dict.items():
        af_ds = aff_load3(f"Data/af/{ds}.af")
        aff_remove_short(af_ds, cut=15)
        aff_load_caid_scores(
            af_ds, f"Data/scores/", prd_list=[score_tag],
            merged=False, remove_missing_scores=True
        )

        ds_data = get_tags_list_scores(
            af_ds, tags=[[tag, "0"], [tag, "1"]], score_tag=score_tag
        )

        labels_list.extend([f"{ds}-0", f"{ds}-1"])
        data.extend(ds_data)

    # Generate plot
    #     t_name = f"{_p}results/Tables/Table_x4_Left.tsv"
    violin_h_plot(
        data=data,
        labels=labels_list,
        display_means={"0": "#d62728", "1": "#2ca02c"},   # red / green
        xlabel="IUPred3 Disorder Scores",
        f_name="Data/results/Figure_4/Figure_4A_left_IUPred.png",
        t_name="Data/results/Tables/Table_4_Left.tsv",
        fontsize=24,
    )
