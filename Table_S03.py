#!/usr/bin/env python3
"""
Supplementary Table S3 — Training AUC ranks for CNN models
across DBsh, CAID1uh, and CAID23uh.

Author: Nawar Malhis
The University of British Columbia, 2026
"""

from param import aff_path
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)


def load_training(in_file):
    """Parse training AUC TSV file."""
    dta = {
        'DBsh': {'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []},
        'CAID1uh': {'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []},
        'CAID23uh': {'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []}
    }

    labels_list = ['Rank', 'CAID1uh', 'CAID23uh', 'DBsh']
    mdl = ''

    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) < 2:
                continue

            lst = line.split()

            if line[0] == '#':
                mdl = lst[1]
                continue

            if lst[0] == 'Rank':
                continue

            for ii, lbl in enumerate(labels_list):
                if ii == 0:
                    continue
                if ii >= len(lst) or lst[ii] == 'NAN':
                    continue
                try:
                    dta[lbl][mdl].append(float(lst[ii]))
                except (ValueError, KeyError):
                    continue

    return dta


if __name__ == '__main__':
    _p = 'Data/'
    dta = load_training(f"{_p}Training_AUC/training_auc.tsv")

    models_list = ['CNN_C1u', 'CNN_C23u', 'CNN_DBs', 'CNN_TR08u']

    out_file = 'Data/results/Tables/Table_S3.tsv'

    with open(out_file, 'w') as fout:
        for mdl in models_list:
            print(f"Model: {mdl}", file=fout)

            # Datasets with sufficient data for this model
            active_datasets = [dd for dd in dta if len(dta[dd][mdl]) > 2]

            print(active_datasets, file=fout)

            # Print ranks (starting from rank 20)
            for jj in range(21):
                rank = jj + 20
                print(rank, end='\t', file=fout)
                for dd in active_datasets:
                    print(dta[dd][mdl][jj], end='\t', file=fout)
                print(file=fout)

            print(file=fout)   # separator between models

    print(f"Table saved to: {out_file}")