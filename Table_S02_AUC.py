#!/usr/bin/env python3
"""
Supplementary Table S02 — AUC values for all predictors
across CAID1uh, CAID23uh, and DBsh.

Author: Nawar Malhis
Refined using Grok
The University of British Columbia, 2026
"""

from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *
from annotated_fasta_metrics import *
from annotated_fasta_CAID import *


if __name__ == '__main__':
    prd_list = [
        'AlphaFold-binding', 'ANCHOR-2', 'CNN_C1u', 'CNN_C23u', 'CNN_DBs', 'CNN_TR08u',
        'DeepDISObind-protein', 'DeepDRPBind-protein', 'DisoRDPbind-protein',
        'DRPBind-protein', 'fMoRFpred', 'OPAL', 'MoRFchibi',
        'MoRFchibi-light', 'MoRFchibi-web'
    ]

    files_dict = {
        'CAID1uh': 'binding_protein',
        'CAID23uh': 'binding_protein',
        'DBsh': 'PDB'
    }

    auc_dict = {}

    for fl in files_dict:
        prd_used = dataset_prd_dict[fl]
        tag = files_dict[fl]

        af = aff_load3(in_file=f"{data_path}af/{fl}.af")
        aff_load_caid_scores(
            af,
            scores_path=f"{data_path}scores/",
            prd_list=prd_used,
            merged=False,
            remove_missing_scores=False
        )

        auc_dict[fl], _ = aff_roc(af, tag=tag, prd_list=prd_used, display=False)

    out_file = "Data/results/Tables/Table_S02_AUC.tsv"

    with open(out_file, 'w') as fout:
        print("Predictor\tCAID1uh\tCAID2&3uh\tDBsh", file=fout)
        for prd in prd_list:
            print(f"{prd}", end=':\t', file=fout)
            for fl in files_dict:
                if prd in auc_dict[fl]:
                    print(f"{auc_dict[fl][prd]:.3f}", end='\t', file=fout)
                else:
                    print(end='\t', file=fout)
            print(file=fout)

    print(f"Table saved to: {out_file}")