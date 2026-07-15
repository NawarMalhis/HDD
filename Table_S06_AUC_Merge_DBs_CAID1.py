from annotated_fasta import *
from param import *
from annotated_fasta_metrics import *
from annotated_fasta_CAID import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)


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
                    continue
                for ii in range(sz):
                    if af[fl]['data'][ac]['tags'][trg][ii] == '1':
                        files_dict[fl]['predictors'][prd]['c1'].append(af[fl]['data'][ac]['scores'][prd][ii])
                    elif af[fl]['data'][ac]['tags'][trg][ii] == '0':
                        files_dict[fl]['predictors'][prd]['c0'].append(af[fl]['data'][ac]['scores'][prd][ii])


def generate_output():
    for fl in files_dict:
        auc_dict[fl] = {}
        for prd in prd_used:
            if merged_class == 'c0':
                yy = (['0'] * (len(files_dict[files[0]]['predictors'][prd]['c0']) +
                               len(files_dict[files[1]]['predictors'][prd]['c0'])) +
                      ['1'] * (len(files_dict[fl]['predictors'][prd]['c1'])))
                sc = files_dict[files[0]]['predictors'][prd]['c0'] + files_dict[files[1]]['predictors'][prd]['c0'] + \
                     files_dict[fl]['predictors'][prd]['c1']
            else:
                yy = (['1'] * (len(files_dict[files[0]]['predictors'][prd]['c1']) +
                               len(files_dict[files[1]]['predictors'][prd]['c1'])) +
                      ['0'] * (len(files_dict[fl]['predictors'][prd]['c0'])))
                sc = files_dict[files[0]]['predictors'][prd]['c1'] + files_dict[files[1]]['predictors'][prd]['c1'] + \
                     files_dict[fl]['predictors'][prd]['c0']
            auc_dict[fl][prd] = roc_auc_score(yy, sc)


if __name__ == '__main__':
    af = {}
    sc_yy = {}
    auc_dict = {}
    files = ['CAID1uh', 'DBsh']
    prd_used = ['ANCHOR-2', 'CNN_C23u', 'CNN_TR08u', 'DisoRDPbind-protein', 'fMoRFpred', 'OPAL', 'MoRFchibi',
                'MoRFchibi-light', 'MoRFchibi-web']

    files_dict = {'CAID1uh': {'tag': 'binding_protein', 'predictors': {}},
                  'DBsh': {'tag': 'PDB', 'predictors': {}}
                  }

    load_data()
    extract_classes()
    for merged_class in ['c0', 'c1']:
        generate_output()
        with open(f"Data/results/Tables/Table_S06_Merge_{merged_class}_{files[0]}_{files[1]}.tsv", 'w') as fout:
            print(f"\nClass_merged {merged_class}\t{files[0]}\t{files[1]}\tDelta", file=fout)
            for prd in prd_used:
                auc_f0 = auc_dict[files[0]][prd]
                auc_f1 = auc_dict[files[1]][prd]
                print(f"{prd}\t{auc_f0:0.3f}\t{auc_f1:0.3f}\t{auc_f0 - auc_f1:0.3f}", flush=True, file=fout)

