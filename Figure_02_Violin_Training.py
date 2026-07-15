import matplotlib.pyplot as plt
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


def violin_plot(data, labels, positions=None, title=None, ax = None,
                    fontsize=18, ylabel='', hh_lines=None):
    ax.violinplot(data, points=200, positions=positions, widths=widths, showmeans=True)
    if title:
        ax.set_title(title, y=0.877, fontsize=fontsize-4)

    for ii, vv in enumerate(hh_lines):
        if vv is None:
            continue
        ax.plot([ii + 0.55, ii + 1.45], [vv, vv], color='red', linewidth=2, linestyle='dashed')

    ax.yaxis.grid(True)
    ax.set_xticks([y + 1 for y in range(len(labels))],
                  labels=labels, fontsize=fontsize-3)
    ax.set_ylabel(ylabel, fontsize=fontsize-4)


if __name__ == '__main__':
    _p = 'Data/'
    dta = load_training(f"{_p}Training_AUC/training_auc.tsv")
    models_list = ['CNN_C1u', 'CNN_C23u', 'CNN_DBs', 'CNN_TR08u']
    dta_list = list(dta.keys())
    plt.rcParams.update({'font.size': 14})
    fig, ax = plt.subplots(nrows=3, ncols=1)
    fig.set_size_inches(12, 10)
    fig.subplots_adjust(left=0.078, right=0.97, top=0.98, bottom=0.05, wspace=0.2, hspace=0.05)

    for ii, dd in enumerate(dta):
        data = []
        widths = []
        labels_list = []
        pos_list = []
        used_auc = []
        ax[ii].set_ylim(0.46, 0.85)
        ax[ii].set_xlim(0.5, 4.5)
        pos = 1

        for mdl in models_list:
            if len(dta[dd][mdl]) > 15:
                used_auc.append(dta[dd][mdl][10])
            else:
                used_auc.append(0.0)
            pos_list.append(pos)
            pos += 1
            data.append(dta[dd][mdl])
            widths.append(0.9)
        if dd == 'CAID23uh':
            labels_list = models_list
        violin_plot(data=data, labels=labels_list, fontsize=22, ylabel=f'{dd} AUC', hh_lines=used_auc,
                        positions=pos_list, ax=ax[ii])
    f_name = f"{_p}results/Figure_2/train_violin_Run_M.png"

    plt.savefig(f_name)
    plt.show()