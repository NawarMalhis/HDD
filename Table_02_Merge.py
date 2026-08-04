#!/usr/bin/env python3
"""
Table 2 — Average ΔAUC under different class-merging strategies:
- non: no merging
- c0: merge Class 0 (non-binding) from both datasets
- c1: merge Class 1 (binding) from both datasets

Author: Nawar Malhis
Refined with assistance from Grok
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
    for fl in files2_list:
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
    for fl in files2_list:
        tag = files_dict_config[fl]
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
    for dst in range(len(files2_list)):
        auc_dict[files2_list[dst]] = {}
        for prd in prd_used:
            scores = {"c0": {}, "c1": {}}
            for cc in [0, 1]:
                scores[f"c{cc}"][files2_list[0]] = files_dict[files2_list[0]]["predictors"][prd][f"c{cc}"]
                scores[f"c{cc}"][files2_list[1]] = files_dict[files2_list[1]]["predictors"][prd][f"c{cc}"]

            if merged_class == "c0":
                # Merge Class 0 from both, keep Class 1 from current dataset
                yy = (["0"] * (len(scores["c0"][files2_list[0]]) + len(scores["c0"][files2_list[1]])) +
                      ["1"] * len(scores["c1"][files2_list[dst]]))
                sc = scores["c0"][files2_list[0]] + scores["c0"][files2_list[1]] + scores["c1"][files2_list[dst]]
            elif merged_class == "c1":
                # Merge Class 1 from both, keep Class 0 from current
                yy = (["1"] * (len(scores["c1"][files2_list[0]]) + len(scores["c1"][files2_list[1]])) +
                      ["0"] * len(scores["c0"][files2_list[dst]]))
                sc = scores["c1"][files2_list[0]] + scores["c1"][files2_list[1]] + scores["c0"][files2_list[dst]]
            else:  # 'non' - no merging
                yy = (["1"] * (len(scores["c1"][files2_list[dst]])) +
                      ["0"] * len(scores["c0"][files2_list[dst]]))
                sc = scores["c1"][files2_list[dst]] + scores["c0"][files2_list[dst]]

            auc_dict[files2_list[dst]][prd] = roc_auc_score(yy, sc)
    return auc_dict


if __name__ == "__main__":
    # Configuration
    files = [["CAID1uh", "DBsh"], ["CAID23uh", "DBsh"], ['CAID1uh', 'CAID23uh']]

    files_dict_config = {
        "CAID1uh": "binding_protein",
        "CAID23uh": "binding_protein",
        "DBsh": "PDB",
    }

    prd_included = {"CAID1uh": ["ANCHOR-2", "CNN_C23u", "CNN_DBs", "CNN_TR08u", "DisoRDPbind-protein", "fMoRFpred",
                                "OPAL", "MoRFchibi", "MoRFchibi-light", "MoRFchibi-web"],

                    "CAID23uh": ["AlphaFold-binding", "ANCHOR-2", "CNN_C1u", 'CNN_DBs', "CNN_TR08u",
                                 "DeepDISObind-protein", "DeepDRPBind-protein", "DisoRDPbind-protein",
                                 "DRPBind-protein", "fMoRFpred", "OPAL", "MoRFchibi", "MoRFchibi-light",
                                 "MoRFchibi-web"],

                    "DBsh": ["AlphaFold-binding", "ANCHOR-2", "CNN_C1u", 'CNN_C23u', "CNN_TR08u",
                             "DeepDISObind-protein", "DeepDRPBind-protein", "DisoRDPbind-protein", "DRPBind-protein",
                             "fMoRFpred", "OPAL", "MoRFchibi", "MoRFchibi-light", "MoRFchibi-web"]
                    }

    out_file = "Data/results/Tables/Table_2_merge.tsv"
    with open(out_file, "w", encoding="utf-8") as fout:
        for fl_i in range(len(files)):  # [0, 1, 2]
            files2_list = files[fl_i]
            prd_used = list(set(prd_included[files2_list[0]]) & set(prd_included[files2_list[1]]))
            print(prd_used)
            # Load and extract
            af = load_data()
            files_dict = extract_class_scores(af)

            # Compute for all merging strategies
            results_dict = {"non": {}, "c0": {}, "c1": {}}  # {files2_list[0]:{}, files2_list[1]:{}}
            for merged_class in ["non", "c0", "c1"]:
                results_dict[merged_class] = compute_auc_merged(files_dict, merged_class)

            # Group definitions
            groups = {
                "AlphaFold-binding": "A", "CNN_C1u": "A", "CNN_C23u": "A", "DeepDISObind-protein": "A",
                "DeepDRPBind-protein": "A", "DisoRDPbind-protein": "A", "DRPBind-protein": "A",
                "fMoRFpred": "B", "MoRFchibi": "B", "OPAL": "B", "CNN_TR08u": "B", "CNN_DBs": "B",
                "MoRFchibi-light": "C", "MoRFchibi-web": "C", "ANCHOR-2": " ",
            }

            # Compute average ΔAUC per group
            avg_groups = {m: {g: {"AUC": 0.0, "cnt": 0} for g in " ABC"} for m in results_dict}

            for merged_class in results_dict:
                for prd in prd_used:
                    auc_c = results_dict[merged_class][files2_list[0]][prd]
                    auc_d = results_dict[merged_class][files2_list[1]][prd]
                    delta = auc_c - auc_d
                    gp = groups.get(prd, " ")
                    avg_groups[merged_class][gp]["AUC"] += delta
                    avg_groups[merged_class][gp]["cnt"] += 1

            # Write to table
            print(f"Average ΔAUC ({files2_list[0]} - {files2_list[1]})", file=fout)
            print("merged_class\tANCHOR-2\tGroup A\tGroup B\tGroup C", file=fout)

            for merged_class in ["non", "c0", "c1"]:
                print(merged_class, end="\t", file=fout)
                for gp in [" ", "A", "B", "C"]:
                    stats = avg_groups[merged_class][gp]
                    avg_delta = stats["AUC"] / stats["cnt"] if stats["cnt"] > 0 else 0.0
                    print(f"{avg_delta:.3f}", end="\t", file=fout)
                print(file=fout)

    print(f"Table saved to: {out_file}")