from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *
from annotated_fasta_metrics import *
from annotated_fasta_CAID import *


def fill_class_0():
    for tg0 in roc_data_dict:
        tag = roc_data_dict[tg0]['tag']
        for ac in af['CAID23uh']['data']:
            sz = len(af['CAID23uh']['data'][ac]['seq'])
            for ii in range(sz):
                if af['CAID23uh']['data'][ac]['tags']['IDR-CAID'][ii] == tag and \
                        af['CAID23uh']['data'][ac]['tags']['binding_protein'][ii] != '1':
                    for prd in af['CAID23uh']['data'][ac]['scores']:
                        if prd not in roc_data_dict[tg0]['CAID23uh']:
                            roc_data_dict[tg0]['CAID23uh'][prd] = {'sc': [], 'mask': ''}
                            roc_data_dict[tg0]['DBsh'][prd] = {'sc': [], 'mask': ''}

                        roc_data_dict[tg0]['CAID23uh'][prd]['sc'].append(af['CAID23uh']['data'][ac]['scores'][prd][ii])
                        roc_data_dict[tg0]['DBsh'][prd]['sc'].append(af['CAID23uh']['data'][ac]['scores'][prd][ii])

    for tg0 in roc_data_dict:
        for df in files_dict:
            for prd in roc_data_dict[tg0][df]:
                sz = len(roc_data_dict[tg0][df][prd]['sc'])
                roc_data_dict[tg0][df][prd]['mask'] = '0' * sz


def fill_class_1():
    for fl in files_dict:
        tag = files_dict[fl]
        for ac in af[fl]['data']:
            for ii in range(len(af[fl]['data'][ac]['seq'])):
                if af[fl]['data'][ac]['tags'][tag][ii] == '1':
                    for prd in af[fl]['data'][ac]['scores']:
                        for tg0 in roc_data_dict:
                            if prd not in roc_data_dict[tg0][fl]:
                                continue
                            roc_data_dict[tg0][fl][prd]['sc'].append(af[fl]['data'][ac]['scores'][prd][ii])
    for tg0 in roc_data_dict:
        for fl in files_dict:
            for prd in roc_data_dict[tg0][fl]:
                sz_sc = len(roc_data_dict[tg0][fl][prd]['sc'])
                sz_msk = len(roc_data_dict[tg0][fl][prd]['mask'])
                roc_data_dict[tg0][fl][prd]['mask'] = roc_data_dict[tg0][fl][prd]['mask'] + '1' * (sz_sc - sz_msk)


if __name__ == '__main__':
    prd_used = ['AlphaFold-binding', 'ANCHOR-2', 'CNN_C1u', 'CNN_TR08u', 'DeepDISObind-protein', 'DeepDRPBind-protein',
                'DisoRDPbind-protein', 'DRPBind-protein', 'fMoRFpred', 'MoRFchibi', 'MoRFchibi-light', 'MoRFchibi-web',
                'OPAL']

    files_dict = {'CAID23uh': 'binding_protein', 'DBsh': 'PDB'}
    roc_data_dict = {'IDR': {'tag': '1', 'CAID23uh': {}, 'DBsh': {}},
                     'PDB': {'tag': '0', 'CAID23uh': {}, 'DBsh': {}}}
    af = {}
    for fl in files_dict:
        af[fl] = aff_load3(in_file=f"{data_path}af/{fl}.af")
        aff_load_caid_scores(af[fl], scores_path=f"{data_path}scores/", prd_list=prd_used, merged=False,
                             remove_missing_scores=False)

    fill_class_0()
    fill_class_1()

    with open("Data/results/Tables/Table_S04_AUC_PDB_IDR.tsv", 'w') as fout:
        for tg0 in roc_data_dict:
            for df in files_dict:
                for prd in roc_data_dict[tg0][df]:
                    yy = list(roc_data_dict[tg0][df][prd]['mask'])
                    auc = roc_auc_score(yy, roc_data_dict[tg0][df][prd]['sc'])
                    print(f"{tg0}\t{df}\t{prd}\t{auc:0.3f}", file=fout)

