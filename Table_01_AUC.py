from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *
from annotated_fasta_metrics import *
from annotated_fasta_CAID import *


if __name__ == '__main__':
    prd_dict = {'all': ['AlphaFold-binding', 'ANCHOR-2', 'DeepDISObind-protein', 'DeepDRPBind-protein',
                        'DisoRDPbind-protein', 'DRPBind-protein', 'fMoRFpred', 'MoRFchibi', 'MoRFchibi-light',
                        'MoRFchibi-web', 'OPAL'],
                'sub': ['ANCHOR-2', 'DisoRDPbind-protein', 'fMoRFpred', 'MoRFchibi', 'MoRFchibi-light',
                        'MoRFchibi-web', 'OPAL']}
    prd_list = ['AlphaFold-binding', 'ANCHOR-2', 'DeepDISObind-protein', 'DeepDRPBind-protein', 'DisoRDPbind-protein',
                'DRPBind-protein', 'fMoRFpred', 'OPAL', 'MoRFchibi', 'MoRFchibi-light', 'MoRFchibi-web']
    files_dict = {'CAID1u': 'binding_protein', 'CAID23u': 'binding_protein', 'DBs': 'PDB'}
    auc_dict = {}

    for fl in files_dict:
        all_sub = 'all'
        if fl == 'CAID1u':
            all_sub = 'sub'
        tag = files_dict[fl]
        af = aff_load3(in_file=f"{data_path}af/{fl}.af")
        aff_load_caid_scores(af, scores_path=f"{data_path}scores/", prd_list=prd_dict[all_sub], merged=False,
                             remove_missing_scores=False)
        auc_dict[fl], _ = aff_roc(af, tag=tag, prd_list=prd_dict[all_sub], display=False)

    with open('Data/results/Tables/Table_1_AUC.tsv', 'w') as fout:
        print("Predictor\tCAID1u\tCAID2&3u\tDBs\t(CAID2&3u-CAID1u)\t(DBs-CAID1u)\t(DBs-CAID2&3u)", file=fout)
        for prd in prd_list:
            print(f"{prd}", end=':\t', file=fout)
            for fl in files_dict:
                if prd in auc_dict[fl]:
                    print(f"{auc_dict[fl][prd]:.3f}", end='\t', file=fout)
                else:
                    print(end='\t', file=fout)

            if prd in auc_dict['CAID1u']:
                print(f"{auc_dict['CAID23u'][prd]-auc_dict['CAID1u'][prd]:.3f}", end='\t', file=fout)
                print(f"{auc_dict['DBs'][prd] - auc_dict['CAID1u'][prd]:.3f}", end='\t', file=fout)
            else:
                print(end='\t\t', file=fout)
            print(f"{auc_dict['DBs'][prd] - auc_dict['CAID23u'][prd]:.3f}", file=fout)
