from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *
from annotated_fasta_metrics import *
from annotated_fasta_CAID import *


def load_data():
    for fl in files:
        print('Loading ', fl, flush=True)
        af[fl] = aff_load3(in_file=f"{data_path}af/{fl}.af")
        aff_load_caid_scores(af[fl], scores_path=f"{data_path}scores/", prd_list=prd_used, merged=False,
                             remove_missing_scores=False)

def extract_classes():
    for fl in files:
        for prd in prd_used:
            files_dict[fl]['predictors'][prd] = {'c0': [], 'c1': []}

    for fl in files:
        trg = files_dict[fl]['tag']
        for ac in af[fl]['data']:
            sz = len(af[fl]['data'][ac]['seq'])
            for prd in prd_used:
                if prd not in af[fl]['data'][ac]['scores']:
                    # print(f"{fl}:{ac}, {prd} not found", flush=True)
                    continue
                for ii in range(sz):
                    if af[fl]['data'][ac]['tags'][trg][ii] == '1':
                        files_dict[fl]['predictors'][prd]['c1'].append(af[fl]['data'][ac]['scores'][prd][ii])
                    elif af[fl]['data'][ac]['tags'][trg][ii] == '0':
                        files_dict[fl]['predictors'][prd]['c0'].append(af[fl]['data'][ac]['scores'][prd][ii])


def generate_output():
    auc_dict = {}
    for fl in files_dict:
        auc_dict[fl] = {}
        for prd in prd_used:
            if merged_class == 'c0':
                yy = (['0'] * (len(files_dict[files[0]]['predictors'][prd]['c0']) +
                               len(files_dict[files[1]]['predictors'][prd]['c0'])) +
                      ['1'] * (len(files_dict[fl]['predictors'][prd]['c1'])))
                sc = files_dict[files[0]]['predictors'][prd]['c0'] + files_dict[files[1]]['predictors'][prd]['c0'] + \
                     files_dict[fl]['predictors'][prd]['c1']
            elif merged_class == 'c1':
                yy = (['1'] * (len(files_dict[files[0]]['predictors'][prd]['c1']) +
                               len(files_dict[files[1]]['predictors'][prd]['c1'])) +
                      ['0'] * (len(files_dict[fl]['predictors'][prd]['c0'])))
                sc = files_dict[files[0]]['predictors'][prd]['c1'] + files_dict[files[1]]['predictors'][prd]['c1'] + \
                     files_dict[fl]['predictors'][prd]['c0']
            else:
                yy = (['1'] * len(files_dict[fl]['predictors'][prd]['c1']) +
                      ['0'] * (len(files_dict[fl]['predictors'][prd]['c0'])))
                sc = files_dict[fl]['predictors'][prd]['c1'] + files_dict[fl]['predictors'][prd]['c0']


            auc_dict[fl][prd] = roc_auc_score(yy, sc)
    return auc_dict


if __name__ == '__main__':
    af = {}
    sc_yy = {}
    # auc_dict = {}
    files = ['CAID23uh', 'DBsh']
    prd_used = ['AlphaFold-binding', 'ANCHOR-2', 'CNN_C1u', 'CNN_TR08u', 'DeepDISObind-protein', 'DeepDRPBind-protein',
                'DisoRDPbind-protein', 'DRPBind-protein', 'fMoRFpred', 'OPAL', 'MoRFchibi', 'MoRFchibi-light',
                'MoRFchibi-web']

    files_dict = {  # 'CAID1uh': {'tag': 'binding_protein', 'predictors': {}},
        'CAID23uh': {'tag': 'binding_protein', 'predictors': {}},
        'DBsh': {'tag': 'PDB', 'predictors': {}}
    }



    load_data()
    extract_classes()
    results_dict = {'non': {}, 'c0': {}, 'c1': {}}
    for merged_class in results_dict:
        results_dict[merged_class] = generate_output()

    groups = {'AlphaFold-binding': 'A', 'CNN_C1u': 'A', 'DeepDISObind-protein': 'A', 'DeepDRPBind-protein': 'A',
              'DisoRDPbind-protein': 'A', 'DRPBind-protein': 'A', 'fMoRFpred': 'B', 'MoRFchibi': 'B', 'OPAL': 'B',
              'CNN_TR08u': 'B', 'MoRFchibi-light': 'C', 'MoRFchibi-web': 'C', 'ANCHOR-2': ' '}
    avg_groups = {'non': {'A': {'AUC': 0, 'cnt': 0}, 'B': {'AUC': 0, 'cnt': 0}, 'C': {'AUC': 0, 'cnt': 0},
                          ' ': {'AUC': 0, 'cnt': 0}},
                  'c0': {'A': {'AUC': 0, 'cnt': 0}, 'B': {'AUC': 0, 'cnt': 0}, 'C': {'AUC': 0, 'cnt': 0},
                         ' ': {'AUC': 0, 'cnt': 0}},
                  'c1': {'A': {'AUC': 0, 'cnt': 0}, 'B': {'AUC': 0, 'cnt': 0}, 'C': {'AUC': 0, 'cnt': 0},
                         ' ': {'AUC': 0, 'cnt': 0}}}

    with open('Data/results/Tables/Table_2_merge.tsv', 'w') as fout:
        print("Average delta AUC values (CAID23uh - DBsh)", file=fout)
        print("merged_class\tANCHOR-2\tA\tB\tC", file=fout)
        # for merged_class in results_dict:
        for merged_class in results_dict:
            for prd in prd_used:
                auc_c = results_dict[merged_class]['CAID23uh'][prd]
                auc_d = results_dict[merged_class]['DBsh'][prd]
                gp = groups[prd]
                avg_groups[merged_class][gp]['AUC'] += (auc_c - auc_d)
                avg_groups[merged_class][gp]['cnt'] += 1
        #
        for merged_class in results_dict:
            print(f"{merged_class}", end="\t", file=fout)
            for gp in [' ', 'A', 'B', 'C']:
                auc = avg_groups[merged_class][gp]['AUC'] / avg_groups[merged_class][gp]['cnt']
                print(f"{auc:0.3f}", end="\t", file=fout)
            print(file=fout)
