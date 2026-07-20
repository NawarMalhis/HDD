#!/usr/bin/env python3
"""
Table 3 — AUC values for short (<=70 AA) and long (>70 AA) binding sites
across CAID1uh, CAID23uh, and DBsh (generated separately for each length category).

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

    # Predictor lists per dataset (as used in your main scripts)
    dataset_prd_dict = {
        'CAID1uh': ['ANCHOR-2', 'DisoRDPbind-protein', 'fMoRFpred', 'MoRFchibi', 'MoRFchibi-light',
                    'MoRFchibi-web', 'OPAL', 'CNN_C23u', 'CNN_DBs', 'CNN_TR08u'],
        'CAID23uh': ['AlphaFold-binding', 'ANCHOR-2', 'DeepDISObind-protein',
                     'DeepDRPBind-protein', 'DisoRDPbind-protein', 'DRPBind-protein',
                     'fMoRFpred', 'MoRFchibi', 'MoRFchibi-light', 'MoRFchibi-web',
                     'OPAL', 'CNN_C1u', 'CNN_DBs', 'CNN_TR08u'],
        'DBsh': ['AlphaFold-binding', 'ANCHOR-2', 'DeepDISObind-protein',
                 'DeepDRPBind-protein', 'DisoRDPbind-protein', 'DRPBind-protein',
                 'fMoRFpred', 'MoRFchibi', 'MoRFchibi-light', 'MoRFchibi-web',
                 'OPAL', 'CNN_C23u', 'CNN_C1u', 'CNN_TR08u']
    }

    for sl in ['_short', '_long']:
        auc_dict = {}
        for fl in files_dict:
            prd_used = dataset_prd_dict[fl]
            tag = files_dict[fl]

            af = aff_load3(in_file=f"{data_path}af/{fl}.af")
            aff_load_caid_scores(af, scores_path=f"{data_path}scores/",
                                 prd_list=prd_used, merged=False,
                                 remove_missing_scores=False)

            # Mask to keep only short or long binding sites
            if sl == '_long':
                aff_tag_size(af, tag, sz_range=[71, 10000])
            else:
                aff_tag_size(af, tag, sz_range=[0, 71])

            auc_dict[fl], _ = aff_roc(af, tag=tag, prd_list=prd_used, display=False)

        # Write table
        out_name = f'Data/results/Tables/Table_3{sl}.tsv'
        with open(out_name, 'w') as fout:
            print("Predictor\tGroup\tCAID1uh\tCAID2&3uh\tDBsh", file=fout)
            for prd in prd_list:
                print(f"{prd}", end=':\t', file=fout)
                for fl in files_dict:
                    if prd in auc_dict[fl]:
                        print(f"{auc_dict[fl][prd]:.3f}", end='\t', file=fout)
                    else:
                        print(end='\t', file=fout)
                print(file=fout)

        print(f"Table written: {out_name}")