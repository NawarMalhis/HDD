#!/usr/bin/env python3
"""
Supplementary Table S05 — AUC under class-merging strategies
(CAID23uh vs DBsh) for Class 0 and Class 1 merging.

Author: Nawar Malhis
Refined: Grok
The University of British Columbia, 2026
"""

from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *
from annotated_fasta_metrics import *
from annotated_fasta_CAID import *


def load_data():
    """Load annotated fasta files and predictor scores."""
    for fl in files:
        print(f"Loading {fl} ...", flush=True)
        af[fl] = aff_load3(in_file=f"{data_path}af/{fl}.af")
        aff_load_caid_scores(
            af[fl],
            scores_path=f"{data_path}scores/",
            prd_list=prd_used,
            merged=False,
            remove_missing_scores=False
        )


def extract_classes():
    """Extract Class 0 and Class 1 scores for every predictor."""
    for fl in files:
        for prd in prd_used:
            files_dict[fl]['predictors'][prd] = {'c0': [], 'c1': []}

    for fl in files:
        trg = files_dict[fl]['tag']
        for ac in af[fl]['data']:
            seq_len = len(af[fl]['data'][ac]['seq'])
            for prd in prd_used:
                if prd not in af[fl]['data'][ac]['scores']:
                    continue
                scores = af[fl]['data'][ac]['scores'][prd]
                tags   = af[fl]['data'][ac]['tags'][trg]
                for ii in range(seq_len):
                    if tags[ii] == '1':
                        files_dict[fl]['predictors'][prd]['c1'].append(scores[ii])
                    elif tags[ii] == '0':
                        files_dict[fl]['predictors'][prd]['c0'].append(scores[ii])


def generate_output(merged_class):
    """Compute AUCs under the specified merging strategy."""
    auc_dict = {}
    for fl in files:
        auc_dict[fl] = {}
        for prd in prd_used:
            c0_0 = files_dict[files[0]]['predictors'][prd]['c0']
            c1_0 = files_dict[files[0]]['predictors'][prd]['c1']
            c0_1 = files_dict[files[1]]['predictors'][prd]['c0']
            c1_1 = files_dict[files[1]]['predictors'][prd]['c1']

            if merged_class == 'c0':
                # Merge Class 0 from both datasets; keep Class 1 of current dataset
                yy = ['0'] * (len(c0_0) + len(c0_1)) + ['1'] * len(files_dict[fl]['predictors'][prd]['c1'])
                sc = c0_0 + c0_1 + files_dict[fl]['predictors'][prd]['c1']
            else:  # merged_class == 'c1'
                # Merge Class 1 from both datasets; keep Class 0 of current dataset
                yy = ['1'] * (len(c1_0) + len(c1_1)) + ['0'] * len(files_dict[fl]['predictors'][prd]['c0'])
                sc = c1_0 + c1_1 + files_dict[fl]['predictors'][prd]['c0']

            auc_dict[fl][prd] = roc_auc_score(yy, sc)

    return auc_dict


if __name__ == '__main__':
    files = ['CAID23uh', 'DBsh']

    prd_used = [
        'AlphaFold-binding', 'ANCHOR-2', 'CNN_C1u', 'CNN_TR08u',
        'DeepDISObind-protein', 'DeepDRPBind-protein', 'DisoRDPbind-protein',
        'DRPBind-protein', 'fMoRFpred', 'OPAL', 'MoRFchibi',
        'MoRFchibi-light', 'MoRFchibi-web'
    ]

    files_dict = {
        'CAID23uh': {'tag': 'binding_protein', 'predictors': {}},
        'DBsh':     {'tag': 'PDB',             'predictors': {}}
    }

    af = {}

    # Load data and extract class scores
    load_data()
    extract_classes()

    # Generate tables for both merging strategies
    for merged_class in ['c0', 'c1']:
        auc_dict = generate_output(merged_class)

        out_file = f"Data/results/Tables/Table_S05_Merge_{merged_class}_{files[0]}_{files[1]}.tsv"
        with open(out_file, 'w') as fout:
            print(f"Class_merged {merged_class}\t{files[0]}\t{files[1]}\tDelta", file=fout)
            for prd in prd_used:
                auc_f0 = auc_dict[files[0]][prd]
                auc_f1 = auc_dict[files[1]][prd]
                print(f"{prd}\t{auc_f0:0.3f}\t{auc_f1:0.3f}\t{auc_f0 - auc_f1:0.3f}", file=fout)

        print(f"Table saved to: {out_file}")