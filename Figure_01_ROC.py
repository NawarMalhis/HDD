from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3, aff_save3, aff_save_fasta
from annotated_fasta_metrics import aff_roc, aff_violin_plot
import numpy as np


if __name__ == '__main__':
    _cnn = False
    _p = 'Data/'
    for in_data in ['DBs', 'CAID1u', 'CAID23u']:
        tag = 'binding_protein'
        if 'DBs' in in_data:
            tag = 'PDB'
        for sl in ['']:  # sl stands for '_short', '_long']:
            target_data = f"{in_data}{sl}"
            tools_list = ['ANCHOR-2', 'MoRFchibi', 'OPAL', 'MoRFchibi-light', 'MoRFchibi-web', 'DisoRDPbind-protein',
                          'fMoRFpred']

            af = aff_load3(f"{_p}af/{target_data}.af")  # _UniParc  _extra
            print(len(af['data']))
            if target_data[:7] == 'CAID23u':
                tools_list = tools_list + ['DeepDISObind-protein', 'AlphaFold-binding', 'DeepDRPBind-protein',
                                           'DRPBind-protein']
                if _cnn:
                    tools_list = tools_list + ['CNN_C1u', 'CNN_DBs', 'CNN_TR08u']
                del af['data']['UPI0000136BB7']
                del af['data']['UPI000004023D']
            elif target_data[:6] == 'CAID1u':
                if _cnn:
                    tools_list = tools_list + ['CNN_C23u', 'CNN_DBs', 'CNN_TR08u']
            elif target_data[:3] == 'DBs':
                tools_list = tools_list + ['DeepDISObind-protein', 'AlphaFold-binding', 'DRPBind-protein',
                                           'DeepDRPBind-protein']
                if _cnn:
                    tools_list = tools_list + ['CNN_C23u', 'CNN_C1u', 'CNN_TR08u']

            aff_load_caid_scores(af, f"{_p}scores/", prd_list=tools_list, merged=False,
                                 remove_missing_scores=False)
            figure_file = f"{_p}results/Figure_1/ROC_{target_data}.png"
            if in_data == 'CAID1u':
                figure_file = f"{_p}results/Figure_S1/ROC_{target_data}.png"
            auc_file = None  # f"{_p}results/CNN/AUC_{target_data}.tsv"
            td = target_data
            if in_data[:7] == 'CAID23u':
                if 'h' in in_data:
                    td = f"CAID2&3uh{sl}"
                else:
                    td = f"CAID2&3u{sl}"
            aff_roc(af, tag=tag, prd_list=tools_list, min_auc=0.2, auc_file=auc_file, title=td,
                    figure_file=figure_file, line_format_dict=prd_dict)
