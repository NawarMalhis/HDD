#!/usr/bin/env python3
"""
Supplementary Figure S2 — LIST-S2 conservation scores (Class 0 vs Class 1)
across training and test datasets.

Author: Nawar Malhis
The University of British Columbia, 2026
"""

from param import *
import sys
# Add AFF project path
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3, aff_remove_short
import matplotlib.pyplot as plt


def violin_h_plot(data, labels, display_means=None, title=None,
                  xlabel='', f_name=None, fontsize=24):
    """Horizontal violin plot with mean lines and table output."""
    plt.rcParams.update({'font.size': 14})

    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig.set_size_inches(8, 14)
    fig.subplots_adjust(left=0.23, right=0.95, top=0.95, bottom=0.08)

    parts = ax.violinplot(data, points=200, showmeans=True,
                          side='high', vert=False)

    if display_means is None:
        display_means = {}

    # Color bodies
    for i, pc in enumerate(parts['bodies']):
        color = display_means.get(str(i), 'lightblue')
        pc.set_facecolor(color)

    # Mean dashed lines
    for cc in display_means:
        ii = int(cc)
        if ii < len(data) and data[ii]:
            fmn = sum(data[ii]) / len(data[ii])
            ax.plot([fmn, fmn], [0.2, len(data) + 0.8],
                    color=display_means[cc], linestyle='dashed', linewidth=2)

    ax.yaxis.grid(True)
    ax.set_yticks([y + 1 for y in range(len(labels))])
    ax.set_yticklabels(labels, fontsize=fontsize - 3)
    ax.set_xlabel(xlabel, fontsize=fontsize)

    if title:
        ax.set_title(title, fontsize=fontsize + 1)

    if f_name:
        plt.savefig(f_name, dpi=350, bbox_inches='tight')

    plt.show()


def get_tags_list_scores(af, tags, score_tag):
    """Extract scores for specific tag/class combinations."""
    d_lst = []
    for jj in range(len(tags)):
        sc = []
        for ac in af['data']:
            for ii in range(len(af['data'][ac]['seq'])):
                tag_name = tags[jj][0]
                if af['data'][ac]['tags'][tag_name][ii] == tags[jj][1]:
                    sc.append(af['data'][ac]['scores'][score_tag][ii])
        d_lst.append(sc)
    return d_lst


if __name__ == '__main__':
    _p = './Data/'
    score_tag = 'LIST-S2'

    labels_list = ['PDB', 'IDR']
    d_set_dict = {
        'TR2008u': 'PDB',
        'DBs': 'PDB',
        'CAID23u': 'binding_protein',
        'CAID1u': 'binding_protein'
    }

    # Baseline
    af = aff_load3(in_file=f"{_p}af/DisProt_2025_06_DBs_extra.af")
    aff_remove_short(af, cut=15)
    aff_load_caid_scores(af, f"{_p}scores/", prd_list=[score_tag],
                         merged=False, remove_missing_scores=True)

    pdb_idr_data = get_tags_list_scores(
        af, tags=[['IDR-CAID', '0'], ['IDR-CAID', '1']], score_tag=score_tag
    )
    data = [pdb_idr_data[0], pdb_idr_data[1]]

    # Other datasets
    for ds in d_set_dict:
        tag = d_set_dict[ds]
        af_ds = aff_load3(in_file=f"{_p}af/{ds}.af")
        aff_remove_short(af_ds, cut=15)
        aff_load_caid_scores(af_ds, scores_path=f"{_p}scores/",
                             prd_list=[score_tag], merged=False,
                             remove_missing_scores=True)

        binding_data = get_tags_list_scores(
            af_ds, tags=[[tag, '0'], [tag, '1']], score_tag=score_tag
        )

        labels_list.extend([f"{ds}-0", f"{ds}-1"])
        data.extend(binding_data)

    vf_name = f"{_p}results/Figure_S2/LIST.png"

    violin_h_plot(
        data=data,
        labels=labels_list,
        display_means={'0': 'red', '1': 'green'},
        xlabel=f"{score_tag} Conservation",
        title="LIST-S2 Conservation Scores",
        f_name=vf_name,
        fontsize=24
    )