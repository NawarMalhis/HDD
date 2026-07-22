#!/usr/bin/env python3
"""
Generate Table 1 — AUC values for IDR-protein binding site predictors
across CAID1u, CAID23u (CAID2&3u), and DBs datasets.

Author: Nawar Malhis
Refined with assistance from Grok
The University of British Columbia, 2026
"""

from param import *
import sys
import numpy as np  # for nan check

if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import aff_load3
from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta_metrics import aff_roc


def main():
    # Configuration
    prd_list = [
        "AlphaFold-binding", "ANCHOR-2", "DeepDISObind-protein",
        "DeepDRPBind-protein", "DisoRDPbind-protein", "DRPBind-protein",
        "fMoRFpred", "OPAL", "MoRFchibi", "MoRFchibi-light", "MoRFchibi-web"
    ]

    dataset_config = {
        "CAID1u": {"tag": "binding_protein", "prd_subset": "sub"},
        "CAID23u": {"tag": "binding_protein", "prd_subset": "all"},
        "DBs": {"tag": "PDB", "prd_subset": "all"},
    }

    prd_dict = {
        "all": [
            "AlphaFold-binding", "ANCHOR-2", "DeepDISObind-protein",
            "DeepDRPBind-protein", "DisoRDPbind-protein", "DRPBind-protein",
            "fMoRFpred", "MoRFchibi", "MoRFchibi-light", "MoRFchibi-web", "OPAL"
        ],
        "sub": ["ANCHOR-2", "DisoRDPbind-protein", "fMoRFpred",
                "MoRFchibi", "MoRFchibi-light", "MoRFchibi-web", "OPAL"],
    }

    auc_dict = {}

    for fl, cfg in dataset_config.items():
        tag = cfg["tag"]
        subset = cfg["prd_subset"]

        af = aff_load3(in_file=f"{data_path}af/{fl}.af")
        aff_load_caid_scores(
            af,
            scores_path=f"Data/scores/",
            prd_list=prd_dict[subset],
            merged=False,
            remove_missing_scores=False,
        )

        auc_dict[fl], _ = aff_roc(
            af, tag=tag, prd_list=prd_dict[subset], display=False
        )

    # Write TSV table
    out_file = f"Data/results/Tables/Table_1_AUC.tsv"

    with open(out_file, "w", encoding="utf-8") as fout:
        header = (
            "Predictor\tCAID1u\tCAID2&3u\tDBs\t"
            "(CAID2&3u - CAID1u)\t(DBs - CAID1u)\t(DBs - CAID2&3u)"
        )
        print(header, file=fout)

        for prd in prd_list:
            row = [prd]

            # AUC values
            for fl in ["CAID1u", "CAID23u", "DBs"]:
                auc = auc_dict[fl].get(prd, float("nan"))
                row.append(f"{auc:.3f}" if not np.isnan(auc) else "")

            # Differences
            if prd in auc_dict["CAID1u"]:
                diff23_1 = auc_dict["CAID23u"].get(prd, 0) - auc_dict["CAID1u"].get(prd, 0)
                diff_db_1 = auc_dict["DBs"].get(prd, 0) - auc_dict["CAID1u"].get(prd, 0)
                diff_db_23 = auc_dict["DBs"].get(prd, 0) - auc_dict["CAID23u"].get(prd, 0)
                row.extend([f"{diff23_1:.3f}", f"{diff_db_1:.3f}", f"{diff_db_23:.3f}"])
            else:
                row.extend(["", "", ""])

            print("\t".join(row), file=fout)

    print(f"Table saved to: {out_file}")


if __name__ == "__main__":
    main()