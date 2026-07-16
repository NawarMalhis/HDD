from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3, aff_remove_short, aff_save3
import matplotlib.pyplot as plt
import matplotlib


def violin_h_plot(data, labels, positions=None, showmeans=True, title=None,
                  display_means=None, fontsize=18, xlabel='', hh_lines=None, f_name=None):
    tbl_name = f"{f_name[0:-4]}_tableLIST.tsv"
    with open(tbl_name, 'w') as fout:
        for ii in range(len(labels)):
            print(f"{labels[ii]}\t{sum(data[ii]) / len(data[ii]):0.3f}", file=fout)
            print(f"{labels[ii]}\t{sum(data[ii]) / len(data[ii]):0.3f}")
    plt.rcParams.update({'font.size': 14})
    if display_means is None:
        display_means = {}
    if hh_lines is None:
        hh_lines = []
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig.set_size_inches(8, 14)
    fig.subplots_adjust(left=0.23, right=0.95, top=0.95, bottom=0.08)
    parts = ax.violinplot(data, points=200, positions=positions, showmeans=showmeans, side='high', vert=False)
    #
    for ii, pc in enumerate(parts['bodies']):
        pc.set_facecolor(display_means[f"{ii}"])
        if ii == 1:
            break

    if title:
        ax.set_title(title, fontsize=fontsize)

    for cc in display_means:
        ii = int(cc)
        fmn = sum(data[ii]) / len(data[ii])
        ax.plot([fmn, fmn], [0.2, len(data)+0.8], color=display_means[cc], linestyle='dashed')
        # print(f"{cc}\t{fmn}")

    for ii, vv in enumerate(hh_lines):
        ax.plot([vv, vv], [ii + 0.8, ii + 1.2], color='red', linestyle='dashed')

    ax.yaxis.grid(True)
    ax.set_yticks([y + 1 for y in range(len(labels))],
                  labels=labels)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    # ax.set_ylabel(ylabel, fontsize=fontsize)
    if f_name:
        plt.savefig(f_name)
    plt.show()

def get_tags_list_scores(af, tags, score_tag):
    d_lst = []
    for jj in range(len(tags)):
        sc = []
        for ac in af['data']:
            # print(f"{ac}\t{len(af['data'][ac]['seq']):,}", flush=True)
            for ii in range(len(af['data'][ac]['seq'])):
                tag = tags[jj][0]
                if af['data'][ac]['tags'][tag][ii] == tags[jj][1]:
                    sc.append(af['data'][ac]['scores'][score_tag][ii])
        d_lst.append(sc)
    return d_lst


if __name__ == '__main__':
    # print(matplotlib.__version__)
    _p = './Data/'
    score_tag = 'LIST-S2'  # 'IUPred3', 'LIST-S2'
    labels_list = ['PDB', 'IDR']
    d_set_dict = {'TR2008u': 'PDB', 'DBs': 'PDB', 'CAID23u': 'binding_protein', 'CAID1u': 'binding_protein'}
    af_dict = {}
    binding_data_dict = {}
    af = aff_load3(f"{_p}af/DisProt_2025_06_DBs_extra.af")
    aff_remove_short(af, cut=15)  # remove short sequences with len(seq) < 15 AA
    aff_load_caid_scores(af, f"{_p}scores/", prd_list=[score_tag], merged=False,
                         remove_missing_scores=True)
    # print(len(af['data']))
    pdb_idr_data = get_tags_list_scores(af=af, tags=[['IDR-CAID', '0'], ['IDR-CAID', '1']], score_tag=score_tag)
    data = [pdb_idr_data[0], pdb_idr_data[1]]
    for ds in d_set_dict:
        tag = d_set_dict[ds]
        af_dict[ds] = aff_load3(f"{_p}af/{ds}.af")
        aff_remove_short(af=af_dict[ds], cut=15)  # remove short sequences with len(seq) < 15 AA
        # print(ds, tag, score_tag, len(af_dict[ds]['data']), list(af_dict[ds]['data'].keys()), flush=True)
        aff_load_caid_scores(af=af_dict[ds], scores_path=f"{_p}scores/", prd_list=[score_tag], merged=False,
                             remove_missing_scores=True)
        binding_data_dict[ds] = get_tags_list_scores(af_dict[ds], tags=[[tag, '0'], [tag, '1']], score_tag=score_tag)
        labels_list.append(f"{ds}-0")
        data.append(binding_data_dict[ds][0])
        labels_list.append(f"{ds}-1")
        data.append(binding_data_dict[ds][1])

    vf_name = f"{_p}results/Figure_S2/LIST.png"
    violin_h_plot(data=data, labels=labels_list, display_means={'0': 'red', '1': 'green'},
                  xlabel=f"{score_tag} Conservation", f_name=vf_name, fontsize=24)
