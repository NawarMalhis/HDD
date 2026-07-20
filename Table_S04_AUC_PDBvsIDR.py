#!/usr/bin/env python3
"""
Supplementary Table S04 — AUC values when treating PDB-derived vs IDR-derived
non-binding residues (Class 0) separately against binding sites (Class 1).

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


def fill_class_0():
    """Collect Class 0 (non-binding) scores from CAID23uh for both IDR and PDB backgrounds."""
    for tg0 in roc_data_dict:
        tag_val = roc_data_dict[tg0]['tag']

        for ac in af['CAID23uh']['data']:
            seq_len = len(af['CAID23uh']['data'][ac]['seq'])
            idr_tags = af['CAID23uh']['data'][ac]['tags']['IDR-CAID']
            bp_tags  = af['CAID23uh']['data'][ac]['tags']['binding_protein']

            for ii in range(seq_len):
                # Keep only non-binding residues that match the current background (IDR or PDB)
                if idr_tags[ii] == tag_val and bp_tags[ii] != '1':
                    for prd in af['CAID23uh']['data'][ac]['scores']:
                        if prd not in roc_data_dict[tg0]['CAID23uh']:
                            roc_data_dict[tg0]['CAID23uh'][prd] = {'sc': [], 'mask': ''}
                            roc_data_dict[tg0]['DBsh'][prd]     = {'sc': [], 'mask': ''}

                        score = af['CAID23uh']['data'][ac]['scores'][prd][ii]
                        roc_data_dict[tg0]['CAID23uh'][prd]['sc'].append(score)
                        roc_data_dict[tg0]['DBsh'][prd]['sc'].append(score)

    # Initialize masks for Class 0
    for tg0 in roc_data_dict:
        for df in files_dict:
            for prd in roc_data_dict[tg0][df]:
                sz = len(roc_data_dict[tg0][df][prd]['sc'])
                roc_data_dict[tg0][df][prd]['mask'] = '0' * sz


def fill_class_1():
    """Append Class 1 (binding) scores from both datasets and update masks."""
    for fl in files_dict:
        tag = files_dict[fl]

        for ac in af[fl]['data']:
            seq_len = len(af[fl]['data'][ac]['seq'])
            for ii in range(seq_len):
                if af[fl]['data'][ac]['tags'][tag][ii] == '1':
                    for prd in af[fl]['data'][ac]['scores']:
                        for tg0 in roc_data_dict:
                            if prd not in roc_data_dict[tg0][fl]:
                                continue
                            score = af[fl]['data'][ac]['scores'][prd][ii]
                            roc_data_dict[tg0][fl][prd]['sc'].append(score)

    # Extend masks with Class 1 labels
    for tg0 in roc_data_dict:
        for fl in files_dict:
            for prd in roc_data_dict[tg0][fl]:
                sc_len  = len(roc_data_dict[tg0][fl][prd]['sc'])
                mask_len = len(roc_data_dict[tg0][fl][prd]['mask'])
                roc_data_dict[tg0][fl][prd]['mask'] += '1' * (sc_len - mask_len)


if __name__ == '__main__':
    prd_used = [
        'AlphaFold-binding', 'ANCHOR-2', 'CNN_C1u', 'CNN_TR08u',
        'DeepDISObind-protein', 'DeepDRPBind-protein', 'DisoRDPbind-protein',
        'DRPBind-protein', 'fMoRFpred', 'MoRFchibi', 'MoRFchibi-light',
        'MoRFchibi-web', 'OPAL'
    ]

    files_dict = {
        'CAID23uh': 'binding_protein',
        'DBsh': 'PDB'
    }

    roc_data_dict = {
        'IDR': {'tag': '1', 'CAID23uh': {}, 'DBsh': {}},
        'PDB': {'tag': '0', 'CAID23uh': {}, 'DBsh': {}}
    }

    # Load data
    af = {}
    for fl in files_dict:
        af[fl] = aff_load3(in_file=f"{data_path}af/{fl}.af")
        aff_load_caid_scores(
            af[fl],
            scores_path=f"{data_path}scores/",
            prd_list=prd_used,
            merged=False,
            remove_missing_scores=False
        )

    # Build score lists
    fill_class_0()
    fill_class_1()

    # Write results
    out_file = "Data/results/Tables/Table_S04_AUC_PDB_IDR.tsv"

    with open(out_file, 'w') as fout:
        for tg0 in roc_data_dict:
            for df in files_dict:
                for prd in roc_data_dict[tg0][df]:
                    yy = list(roc_data_dict[tg0][df][prd]['mask'])
                    sc = roc_data_dict[tg0][df][prd]['sc']
                    if len(yy) == len(sc) and len(yy) > 0:
                        auc = roc_auc_score(yy, sc)
                        print(f"{tg0}\t{df}\t{prd}\t{auc:0.3f}", file=fout)

    print(f"Table saved to: {out_file}")