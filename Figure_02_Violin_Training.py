#!/usr/bin/env python3
"""
Generate violin plots comparing training AUC distributions for CNN models
across DBsh, CAID1uh, and CAID23uh datasets.

Author: Nawar Malhis
The University of British Columbia, 2026
"""

import matplotlib.pyplot as plt
from pathlib import Path
from param import *
import sys

# Add AFF path
if aff_path not in globals() or aff_path not in sys.path:  # type: ignore[name-defined]
    sys.path.append(str(aff_path))  # type: ignore[name-defined]


def load_training_auc(in_file: str | Path) -> dict:
    """
    Parse training AUC TSV file into nested dictionary.
    Expected format: header lines with model names, then rank + AUC columns.
    """
    data = {
        "DBsh": {"CNN_C1u": [], "CNN_C23u": [], "CNN_DBs": [], "CNN_TR08u": []},
        "CAID1uh": {"CNN_C1u": [], "CNN_C23u": [], "CNN_DBs": [], "CNN_TR08u": []},
        "CAID23uh": {"CNN_C1u": [], "CNN_C23u": [], "CNN_DBs": [], "CNN_TR08u": []},
    }

    labels_list = ["Rank", "CAID1uh", "CAID23uh", "DBsh"]

    with open(in_file, encoding="utf-8") as fin:
        for line in fin:
            line = line.strip()
            if len(line) < 2:
                continue

            lst = line.split()

            if line.startswith("#"):
                # e.g., # CNN_C1u
                if len(lst) > 1:
                    model = lst[1]
                continue

            if lst[0] == "Rank":
                continue

            # Parse AUC values for each dataset
            for col_idx, dataset in enumerate(labels_list):
                if col_idx == 0:  # Skip Rank column
                    continue
                if col_idx >= len(lst) or lst[col_idx] == "NAN":
                    continue
                try:
                    data[dataset][model].append(float(lst[col_idx]))
                except (ValueError, NameError, KeyError):
                    continue  # Skip malformed entries

    return data


def violin_plot(
    data: list,
    labels: list[str],
    positions: list | None = None,
    title: str | None = None,
    ax: plt.Axes | None = None,
    fontsize: int = 18,
    ylabel: str = "",
    hh_lines: list | None = None,
    widths: list | None = None,
) -> None:
    """Draw a single violin plot with optional horizontal reference lines."""
    if ax is None:
        ax = plt.gca()

    if widths is None:
        widths = [0.9] * len(data)

    ax.violinplot(
        data,
        points=200,
        positions=positions,
        widths=widths,
        showmeans=True,
    )

    if title:
        ax.set_title(title, y=0.88, fontsize=fontsize - 2)

    # Horizontal reference lines (e.g., median or specific rank AUC)
    if hh_lines:
        for i, val in enumerate(hh_lines):
            if val is None or val <= 0:
                continue
            ax.plot(
                [i + 0.55, i + 1.45],
                [val, val],
                color="red",
                linewidth=2,
                linestyle="--",
            )

    ax.yaxis.grid(True, linestyle="--", alpha=0.7)
    ax.set_xticks([y + 1 for y in range(len(labels))])
    ax.set_xticklabels(labels, fontsize=fontsize - 3)
    ax.set_ylabel(ylabel, fontsize=fontsize - 2)
    ax.set_ylim(0.46, 0.85)
    ax.set_xlim(0.5, 4.5)


if __name__ == "__main__":
    BASE_DIR = Path("Data")
    tsv_path = BASE_DIR / "Training_AUC" / "training_auc.tsv"

    dta = load_training_auc(tsv_path)

    models = ["CNN_C1u", "CNN_C23u", "CNN_DBs", "CNN_TR08u"]
    datasets = list(dta.keys())

    plt.rcParams.update({"font.size": 14})
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 10))
    fig.subplots_adjust(left=0.08, right=0.97, top=0.96, bottom=0.08, hspace=0.12)

    for i, dataset in enumerate(datasets):
        ax = axes[i]
        plot_data = []
        positions = []
        widths = []
        ref_lines = []  # AUC values for red dashed lines

        pos = 1
        for model in models:
            aucs = dta[dataset][model]
            plot_data.append(aucs)

            # Use 11th value (0-based index 10) as reference if enough data
            ref = aucs[10] if len(aucs) > 15 else 0.0
            ref_lines.append(ref)

            positions.append(pos)
            widths.append(0.9)
            pos += 1

        # Only show model names on bottom panel
        x_labels = models if dataset == "CAID23uh" else []

        violin_plot(
            data=plot_data,
            labels=x_labels,
            positions=positions,
            ylabel=f"{dataset} AUC",
            hh_lines=ref_lines,
            fontsize=22,
            ax=ax,
        )

    output_dir = BASE_DIR / "results" / "Figure_2"
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / "train_violin_Run_M.png"

    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"Figure saved to: {out_file}")