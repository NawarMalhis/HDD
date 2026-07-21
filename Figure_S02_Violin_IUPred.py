#!/usr/bin/env python3
"""
Supplementary Figure S2 — IUPred3 disorder scores (Class 0 vs Class 1)

Author: Nawar Malhis
Refined using Grok
The University of British Columbia, 2026
"""

import os
os.environ["QT_LOGGING_RULES"] = "qt.qpa.wayland.textinput=false"
import sys
from param import *
# Add AFF project path
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3, aff_remove_short
import matplotlib.pyplot as plt


def get_tags_list_scores(af, tags, score_tag):
    """Extract scores for specific tag/class combinations."""
    d_lst = []
    for jj in range(len(tags)):
        sc = []
        for ac in af['data']:
            for ii in range(len(af['data'][ac]['seq'])):
                tag_name = tags[jj][0]
                if af['data'][ac]['tags'][tag_name][ii] == tags[jj][1]:
                    sc.append(af['data'][ac]['scores'][score_tag][ii])
        d_lst.append(sc)
    return d_lst


def violin_h_plot(
    data: list,
    labels: list,
    display_means: dict | None = None,
    title: str | None = None,
    xlabel: str = "",
    f_name: str | None = None,
    fontsize: int = 24,
):
    """Horizontal violin plot with mean lines."""
    if display_means is None:
        display_means = {}

    # Write means to table
    if f_name:
        tbl_name = f"Data/results/Tables/Table_4_IUPred.tsv"
        with open(tbl_name, "w", encoding="utf-8") as fout:
            for i, lbl in enumerate(labels):
                mean_val = sum(data[i]) / len(data[i]) if data[i] else 0.0
                fout.write(f"{lbl}\t{mean_val:.3f}\n")
                print(f"{lbl}\t{mean_val:.3f}")

    plt.rcParams.update({"font.size": 14})
    fig, ax = plt.subplots(figsize=(8.5, 14))
    fig.subplots_adjust(left=0.26, right=0.96, top=0.96, bottom=0.08)

    positions = list(range(1, len(data) + 1))
    ax.set_ylim(bottom=0.5, top=10.8)
    parts = ax.violinplot(
        data,
        positions=positions,
        points=200,
        showmeans=True,
        vert=False,
        side="high",
    )

    # Color violins
    for i, pc in enumerate(parts["bodies"]):
        color = display_means.get(str(i), "lightblue")
        pc.set_facecolor(color)
        pc.set_edgecolor("black")
        pc.set_alpha(0.85)

    # Mean dashed lines
    for i_str, color in display_means.items():
        i = int(i_str)
        if i < len(data) and data[i]:
            mean_val = sum(data[i]) / len(data[i])
            ax.plot([mean_val, mean_val], [0.2, len(data) + 0.8],
                    color=color, linestyle="--", linewidth=2)

    ax.yaxis.grid(True, linestyle="--", alpha=0.7)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=fontsize - 3)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    if title:
        ax.set_title(title, fontsize=fontsize + 1, pad=20)

    if f_name:
        plt.savefig(f_name, dpi=350, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    score_tag = "IUPred3"   # change to 'LIST-S2' if needed

    af = aff_load3(f"Data/af/DisProt_2025_06_DBs_extra.af")
    aff_remove_short(af, cut=15)
    aff_load_caid_scores(
        af, scores_path=f"Data/scores/", prd_list=[score_tag],
        merged=False, remove_missing_scores=True
    )

    # Base reference data [PDB, IDR]
    labels_list = ["PDB", "IDR"]
    data = get_tags_list_scores(
        af, tags=[["IDR-CAID", "0"], ["IDR-CAID", "1"]], score_tag=score_tag
    )

    d_set_dict = {
        "TR2008u": "PDB",
        "DBs": "PDB",
        "CAID23u": "binding_protein",
        "CAID1u": "binding_protein",
    }

    for ds, tag in d_set_dict.items():
        af_ds = aff_load3(f"Data/af/{ds}.af")
        aff_remove_short(af_ds, cut=15)
        aff_load_caid_scores(
            af_ds, scores_path=f"Data/scores/", prd_list=[score_tag],
            merged=False, remove_missing_scores=True
        )

        ds_data = get_tags_list_scores(
            af_ds, tags=[[tag, "0"], [tag, "1"]], score_tag=score_tag
        )

        labels_list.extend([f"{ds}-0", f"{ds}-1"])
        data.extend(ds_data)

    # Generate plot
    vf_name = f"Data/results/Figure_S2/IUPred3.png"
    violin_h_plot(
        data=data,
        labels=labels_list,
        display_means={"0": "#d62728", "1": "#2ca02c"},  # red / green
        xlabel="IUPred3 Disorder Scores",
        # title="IUPred3 Disorder Scores",
        f_name=vf_name,
        fontsize=24,
    )
