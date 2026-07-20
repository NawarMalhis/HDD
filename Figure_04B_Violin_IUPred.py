#!/usr/bin/env python3
"""
Figure 4B (right) — Horizontal violin plot of IUPred3 scores
for short vs long binding sites (Class 1) across test datasets.

Author: Nawar Malhis
The University of British Columbia, 2026
"""

from param import *
import sys
# Add AFF project path
if aff_path not in sys.path:
    sys.path.append(aff_path)

import matplotlib.pyplot as plt
from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3, aff_remove_short, aff_tag_size


def get_tags_list_scores(af: dict, tags: list, score_tag: str) -> list:
    """Extract scores for given tag/class combinations."""
    data = []
    for tag_name, class_val in tags:
        sc = []
        for entry in af["data"].values():
            tags_dict = entry.get("tags", {})
            scores_dict = entry.get("scores", {})
            if tag_name not in tags_dict or score_tag not in scores_dict:
                continue
            tag_array = tags_dict[tag_name]
            score_array = scores_dict[score_tag]
            for i, t_val in enumerate(tag_array):
                if t_val == class_val:
                    sc.append(score_array[i])
        data.append(sc)
    return data


def violin_h_plot(
    data: list,
    labels: list,
    display_means: dict | None = None,
    title: str | None = None,
    xlabel: str = "",
    f_name: str | None = None,
    t_name: str | None = None,
    fontsize: int = 22,
):
    """Horizontal violin plot with mean lines and table output."""
    if display_means is None:
        display_means = {}

    # Write mean values to TSV
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
        vert=False,
        side="high",
    )

    # Color violin bodies
    for i, pc in enumerate(parts["bodies"]):
        color = display_means.get(str(i), "lightblue")
        pc.set_facecolor(color)
        pc.set_edgecolor("black")
        pc.set_alpha(0.85)

    # Dashed mean lines
    for i_str, color in display_means.items():
        i = int(i_str)
        if i < len(data) and data[i]:
            mean_val = sum(data[i]) / len(data[i])
            ax.plot([mean_val, mean_val], [0.2, len(data) + 0.8],
                    color=color, linestyle="--", linewidth=2.2)

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
    # base_path = Path("Data")
    score_tag = "IUPred3"   # or 'LIST-S2'
    # output_dir = base_path / "results" / "Figure_4"

    # Baseline (DisProt extra)
    af = aff_load3("Data/af/DisProt_2025_06_DBs_extra.af")
    aff_remove_short(af, cut=15)
    aff_load_caid_scores(
        af, scores_path=f"Data/scores/", prd_list=[score_tag],
        merged=False, remove_missing_scores=True
    )

    labels_list = ["PDB", "IDR"]
    data = get_tags_list_scores(
        af, tags=[["IDR-CAID", "0"], ["IDR-CAID", "1"]], score_tag=score_tag
    )

    # Test datasets: short + long binding sites (Class 1 only)
    d_set_dict = {
        "DBsh": "PDB",
        "CAID23uh": "binding_protein",
        "CAID1uh": "binding_protein",
    }
    sl_dict = {"_short": [0, 71], "_long": [71, 10000]}

    for ds, tag in d_set_dict.items():
        af_ds = aff_load3(f"Data/af/{ds}.af")
        aff_remove_short(af_ds, cut=15)

        for sl_label, sz_range in sl_dict.items():
            af_temp = aff_load3(f"Data/af/{ds}.af")  # fresh copy
            aff_remove_short(af_temp, cut=15)
            aff_tag_size(af_temp, tag=tag, sz_range=sz_range)

            aff_load_caid_scores(
                af_temp, scores_path=f"Data/scores/", prd_list=[score_tag],
                merged=False, remove_missing_scores=True
            )

            # Only Class 1 (binding)
            ds_data = get_tags_list_scores(
                af_temp, tags=[[tag, "1"]], score_tag=score_tag
            )

            labels_list.append(f"{ds}{sl_label}")
            data.extend(ds_data)

    print(f"Total violins: {len(data)} | Labels: {len(labels_list)}")

    # Generate plot
    violin_h_plot(
        data=data,
        labels=labels_list,
        display_means={"0": "#d62728", "1": "#2ca02c"},   # red for baseline, green for binding
        xlabel="IUPred3 Disorder Scores",
        f_name="Data/results/Figure_4/Figure_4B_right_IUPred.png",
        t_name="Data/results/Tables/Table_4_Right.tsv",
        fontsize=24,
    )