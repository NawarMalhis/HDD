#!/usr/bin/env python3
"""
Supplementary Table S3 — Training AUC ranks for CNN models
across DBsh, CAID1uh, and CAID23uh.

Author: Nawar Malhis
Refined using Grok
The University of British Columbia, 2026
"""

from param import aff_path
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)


def load_training(in_file):
    """Parse the training AUC TSV file into a nested dictionary."""
    dta = {
        'DBsh':    {'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []},
        'CAID1uh': {'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []},
        'CAID23uh':{'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []}
    }

    labels_list = ['Rank', 'CAID1uh', 'CAID23uh', 'DBsh']
    mdl = ''

    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) < 2:
                continue

            lst = line.split()

            # Model header line, e.g. "# CNN_C1u"
            if line.startswith('#'):
                if len(lst) > 1:
                    mdl = lst[1]
                continue

            # Skip header row
            if lst[0] == 'Rank':
                continue

            # Extract AUC values for each dataset
            for ii, lbl in enumerate(labels_list):
                if ii == 0:                    # skip Rank column
                    continue
                if ii >= len(lst) or lst[ii] == 'NAN':
                    continue
                try:
                    dta[lbl][mdl].append(float(lst[ii]))
                except (ValueError, KeyError):
                    continue

    return dta


if __name__ == '__main__':
    input_file = 'Data/Training_AUC/training_auc.tsv'
    out_file   = 'Data/results/Tables/Table_S3.tsv'

    dta = load_training(input_file)

    models_list = ['CNN_C1u', 'CNN_C23u', 'CNN_DBs', 'CNN_TR08u']

    with open(out_file, 'w') as fout:
        for mdl in models_list:
            print(f"Model: {mdl}", file=fout)

            # Only include datasets that have enough values for this model
            active_datasets = [dd for dd in dta if len(dta[dd][mdl]) > 2]
            print('\t'.join(active_datasets), file=fout)

            # Print ranks 20–40
            for jj in range(21):
                rank = jj + 20
                print(rank, end='\t', file=fout)
                for dd in active_datasets:
                    print(f"{dta[dd][mdl][jj]:.3f}", end='\t', file=fout)
                print(file=fout)

            print(file=fout)   # blank line between models

    print(f"Table saved to: {out_file}")