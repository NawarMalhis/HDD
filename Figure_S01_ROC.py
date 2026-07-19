#!/usr/bin/env python3
"""
Supplementary Figure S1 — ROC curves for DBsh, CAID1uh, CAID23uh

Author: Nawar Malhis
The University of British Columbia, 2026
"""

from pathlib import Path
import sys

# Add AFF project path
from param import *
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3
from annotated_fasta_metrics import aff_roc


def get_tools_list(dataset: str, use_cnn: bool = True) -> list[str]:
    """Return predictor list for a given dataset."""
    base = [
        "ANCHOR-2", "MoRFchibi", "OPAL", "MoRFchibi-light",
        "MoRFchibi-web", "DisoRDPbind-protein", "fMoRFpred"
    ]

    extra = {
        "DBsh": ["DeepDISObind-protein", "AlphaFold-binding",
                 "DRPBind-protein", "DeepDRPBind-protein"],
        "CAID1uh": [],
        "CAID23uh": ["DeepDISObind-protein", "AlphaFold-binding",
                     "DeepDRPBind-protein", "DRPBind-protein"],
    }

    cnn_models = {
        "DBsh": ["CNN_C23u", "CNN_C1u", "CNN_TR08u"],
        "CAID1uh": ["CNN_C23u", "CNN_DBs", "CNN_TR08u"],
        "CAID23uh": ["CNN_C1u", "CNN_DBs", "CNN_TR08u"],
    }

    tools = base + extra.get(dataset, [])
    if use_cnn:
        tools += cnn_models.get(dataset, [])

    return tools


if __name__ == "__main__":
    use_cnn = True
    base_dir = Path("Data")
    results_dir = base_dir / "results" / "Figure_S1"

    datasets = ["DBsh", "CAID1uh", "CAID23uh"]

    for in_data in datasets:
        tag = "PDB" if "DBs" in in_data else "binding_protein"

        for sl in [""]:   # extend with '_short', '_long' when needed
            target_data = f"{in_data}{sl}"

            # Load data
            af_path = base_dir / "af" / f"{target_data}.af"
            af = aff_load3(str(af_path))
            print(f"Loaded {len(af['data'])} sequences for {target_data}")

            # Special cleanup for CAID23uh
            if target_data.startswith("CAID23u"):
                af["data"].pop("UPI0000136BB7", None)
                af["data"].pop("UPI000004023D", None)

            # Load scores
            tools_list = get_tools_list(in_data, use_cnn=use_cnn)
            aff_load_caid_scores(
                af,
                f"{base_dir}/scores/",
                prd_list=tools_list,
                merged=False,
                remove_missing_scores=False,
            )

            # Output figure
            results_dir.mkdir(parents=True, exist_ok=True)
            figure_file = f"{results_dir}/ROC_{target_data}.png"

            # Nice title
            title = target_data
            if in_data.startswith("CAID23u"):
                h = "h" if "h" in in_data else ""
                title = f"CAID2&3u{h}{sl}"

            aff_roc(
                af,
                tag=tag,
                prd_list=tools_list,
                min_auc=0.2,
                auc_file=None,   # set path if you want TSV output
                title=title,
                figure_file=figure_file,
                line_format_dict=prd_dict,  # type: ignore[name-defined]
            )