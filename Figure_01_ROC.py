#!/usr/bin/env python3
"""
Generate ROC curves for protein binding site prediction tools across multiple datasets.
Supports DBs, CAID1u, and CAID23u benchmarks.

Author: Nawar Malhis
The University of British Columbia, 2026
"""

from param import *
import sys

# Add AFF project path
if aff_path not in globals() or aff_path not in sys.path:  # type: ignore[name-defined]
    sys.path.append(str(aff_path))  # type: ignore[name-defined]

from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3
from annotated_fasta_metrics import aff_roc

# ====================== CONFIGURATION ======================
RESULTS_DIR = "Data/results"

# Common tools (baseline predictors)
BASE_TOOLS = [
    "ANCHOR-2",
    "MoRFchibi",
    "OPAL",
    "MoRFchibi-light",
    "MoRFchibi-web",
    "DisoRDPbind-protein",
    "fMoRFpred",
]

# Dataset-specific additional tools
EXTRA_TOOLS = {
    "DBs": ["DeepDISObind-protein", "AlphaFold-binding", "DRPBind-protein", "DeepDRPBind-protein"],
    "CAID1u": [],
    "CAID23u": ["DeepDISObind-protein", "AlphaFold-binding", "DRPBind-protein", "DeepDRPBind-protein"],
}

# Entries to remove from CAID23u (problematic sequences)
CAID23U_TO_REMOVE = {"UPI0000136BB7", "UPI000004023D"}


def get_tools_list(dataset: str) -> list[str]:
    """Build the list of predictors for a given dataset."""
    tools = BASE_TOOLS.copy()

    # Add dataset-specific tools
    tools.extend(EXTRA_TOOLS.get(dataset, []))

    return tools


def main() -> None:
    for in_data in ["DBs", "CAID1u", "CAID23u"]:
        tag = "PDB" if "DBs" in in_data else "binding_protein"

        # Load annotated fasta
        af_path = f"Data/af/{in_data}.af"
        af = aff_load3(str(af_path))
        print(f"Loaded {len(af['data'])} sequences for {in_data}")

        # Remove problematic entries for CAID23u
        if in_data.startswith("CAID23u"):
            for upi in CAID23U_TO_REMOVE:
                af["data"].pop(upi, None)

        # Load scores
        tools_list = get_tools_list(in_data)
        aff_load_caid_scores(
            af,
            scores_path=f"Data/scores/"  ,
            prd_list=tools_list,
            merged=False,
            remove_missing_scores=False,
        )

        # Output paths and title
        if in_data == "CAID1u":
            fig_dir = "Data/results/Figure_S1"
        else:
            fig_dir = "Data/results/Figure_1"

        figure_file = f"{fig_dir}ROC_{in_data}.png"

        title = in_data

        # Generate ROC plot
        aff_roc(
            af,
            tag=tag,
            prd_list=tools_list,
            min_auc=0.2,
            auc_file=None,  # Set to a path if you want TSV output
            title=title,
            figure_file=str(figure_file),
            line_format_dict=prd_dict,  # type: ignore[name-defined]
        )


if __name__ == "__main__":
    main()