from param import aff_path
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)


def load_training(in_file):
    dta = {'DBsh': {'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []},
           'CAID1uh': {'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []},
           'CAID23uh': {'CNN_C1u': [], 'CNN_C23u': [], 'CNN_DBs': [], 'CNN_TR08u': []}}
    mdl = ''
    labels_list = ['Rank', 'CAID1uh', 'CAID23uh', 'DBsh']
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
                if lst[ii] == 'NAN':
                    continue
                dta[lbl][mdl].append(float(lst[ii]))
    return dta


if __name__ == '__main__':
    _p = 'Data/'
    dta = load_training(f"{_p}Training_AUC/training_auc.tsv")
    models_list = ['CNN_C1u', 'CNN_C23u', 'CNN_DBs', 'CNN_TR08u']
    dta_list = list(dta.keys())
    # print(dta_list)
    with open('Data/results/Tables/Table_S3.tsv', 'w') as fout:
        for mdl in models_list:
            print('Model:', mdl, file=fout)
            lst = []
            for dd in dta:
                if len(dta[dd][mdl]) > 2:
                    lst.append(dd)
            print(lst, file=fout)
            for jj in range(21):
                print(jj + 20, end='\t', file=fout)
                for dd in lst:
                    print(dta[dd][mdl][jj], end='\t', file=fout)
                print(file=fout)
            print(file=fout)
