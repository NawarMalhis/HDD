#!/usr/bin/env python3
"""
Table 2 — Average ΔAUC (CAID23uh - DBsh) under different class-merging strategies:
- non: no merging
- c0: merge Class 0 (non-binding) from both datasets
- c1: merge Class 1 (binding) from both datasets

Author: Nawar Malhis
Refined using Grok
The University of British Columbia, 2026
"""

import sys
from sklearn.metrics import roc_auc_score
from param import *

if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import aff_load3
from annotated_fasta_CAID import aff_load_caid_scores


def load_data() -> dict:
    """Load annotated files and predictor scores."""
    af = {}
    for fl in files:
        print(f"Loading {fl} ...", flush=True)
        af[fl] = aff_load3(in_file=f"Data/af/{fl}.af")
        aff_load_caid_scores(
            af[fl],
            scores_path="Data/scores/",
            prd_list=prd_used,
            merged=False,
            remove_missing_scores=False,
        )
    return af


def extract_class_scores(af_dict: dict) -> dict:
    """Extract predictor scores for Class 0 and Class 1 per dataset."""
    files_dict = {}
    for fl in files:
        tag = files_dict_config[fl]["tag"]
        files_dict[fl] = {"tag": tag, "predictors": {}}

        for prd in prd_used:
            files_dict[fl]["predictors"][prd] = {"c0": [], "c1": []}

        for entry in af_dict[fl]["data"].values():
            seq_len = len(entry["seq"])
            trg_tags = entry["tags"][tag]

            for prd in prd_used:
                if prd not in entry.get("scores", {}):
                    continue
                scores = entry["scores"][prd]
                for i in range(seq_len):
                    if trg_tags[i] == "1":
                        files_dict[fl]["predictors"][prd]["c1"].append(scores[i])
                    elif trg_tags[i] == "0":
                        files_dict[fl]["predictors"][prd]["c0"].append(scores[i])
    return files_dict


def compute_auc_merged(files_dict: dict, merged_class: str) -> dict:
    """Compute AUC for each predictor under a specific merging strategy."""
    auc_dict = {}
    for fl in files:
        auc_dict[fl] = {}
        for prd in prd_used:
            c0_23 = files_dict["CAID23uh"]["predictors"][prd]["c0"]
            c1_23 = files_dict["CAID23uh"]["predictors"][prd]["c1"]
            c0_db = files_dict["DBsh"]["predictors"][prd]["c0"]
            c1_db = files_dict["DBsh"]["predictors"][prd]["c1"]

            if merged_class == "c0":
                # Merge Class 0 from both, keep Class 1 from current dataset
                yy = ["0"] * (len(c0_23) + len(c0_db)) + ["1"] * len(c1_23 if fl == "CAID23uh" else c1_db)
                sc = c0_23 + c0_db + (c1_23 if fl == "CAID23uh" else c1_db)
            elif merged_class == "c1":
                # Merge Class 1 from both, keep Class 0 from current
                yy = ["1"] * (len(c1_23) + len(c1_db)) + ["0"] * len(c0_23 if fl == "CAID23uh" else c0_db)
                sc = c1_23 + c1_db + (c0_23 if fl == "CAID23uh" else c0_db)
            else:  # 'non' - no merging
                yy = ["1"] * len(c1_23 if fl == "CAID23uh" else c1_db) + \
                     ["0"] * len(c0_23 if fl == "CAID23uh" else c0_db)
                sc = (c1_23 if fl == "CAID23uh" else c1_db) + \
                     (c0_23 if fl == "CAID23uh" else c0_db)

            auc_dict[fl][prd] = roc_auc_score(yy, sc)
    return auc_dict


if __name__ == "__main__":
    # Configuration
    files = ["CAID23uh", "DBsh"]

    files_dict_config = {
        "CAID23uh": {"tag": "binding_protein"},
        "DBsh": {"tag": "PDB"},
    }

    prd_used = [
        "AlphaFold-binding", "ANCHOR-2", "CNN_C1u", "CNN_TR08u",
        "DeepDISObind-protein", "DeepDRPBind-protein", "DisoRDPbind-protein",
        "DRPBind-protein", "fMoRFpred", "OPAL", "MoRFchibi",
        "MoRFchibi-light", "MoRFchibi-web"
    ]

    # Load and extract
    af = load_data()
    files_dict = extract_class_scores(af)

    # Compute for all merging strategies
    results_dict = {}
    for merged_class in ["non", "c0", "c1"]:
        results_dict[merged_class] = compute_auc_merged(files_dict, merged_class)

    # Group definitions
    groups = {
        "AlphaFold-binding": "A", "CNN_C1u": "A", "DeepDISObind-protein": "A",
        "DeepDRPBind-protein": "A", "DisoRDPbind-protein": "A", "DRPBind-protein": "A",
        "fMoRFpred": "B", "MoRFchibi": "B", "OPAL": "B", "CNN_TR08u": "B",
        "MoRFchibi-light": "C", "MoRFchibi-web": "C", "ANCHOR-2": " ",
    }

    # Compute average ΔAUC per group
    avg_groups = {m: {g: {"AUC": 0.0, "cnt": 0} for g in " ABC"} for m in results_dict}

    for merged_class in results_dict:
        for prd in prd_used:
            auc_c = results_dict[merged_class]["CAID23uh"][prd]
            auc_d = results_dict[merged_class]["DBsh"][prd]
            delta = auc_c - auc_d
            gp = groups.get(prd, " ")
            avg_groups[merged_class][gp]["AUC"] += delta
            avg_groups[merged_class][gp]["cnt"] += 1

    # Write table
    out_file = "Data/results/Tables/Table_2_merge.tsv"

    with open(out_file, "w", encoding="utf-8") as fout:
        print("Average ΔAUC (CAID23uh - DBsh)", file=fout)
        print("merged_class\tANCHOR-2\tGroup A\tGroup B\tGroup C", file=fout)

        for merged_class in ["non", "c0", "c1"]:
            print(merged_class, end="\t", file=fout)
            for gp in [" ", "A", "B", "C"]:
                stats = avg_groups[merged_class][gp]
                avg_delta = stats["AUC"] / stats["cnt"] if stats["cnt"] > 0 else 0.0
                print(f"{avg_delta:.3f}", end="\t", file=fout)
            print(file=fout)

    print(f"Table saved to: {out_file}")