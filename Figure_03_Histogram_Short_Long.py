from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import aff_load3, aff_tag_size, aff_gen_counts
from annotated_fasta_metrics import aff_roc, aff_violin_plot
import numpy as np


if __name__ == '__main__':
    data_dict = {'DisProt': {'tag': 'binding_protein', 'size': [[0, 31], [31, 71], [71, 100000]]},
                  'CAID1u': {'tag': 'binding_protein', 'size': [[0, 31], [31, 71], [71, 100000]]},
                 'CAID23u': {'tag': 'binding_protein', 'size': [[0, 31], [31, 71], [71, 100000]]},
                 'DBs': {'tag': 'PDB', 'size': [[0, 31], [31, 71], [71, 100000]]},
                 'CAID1uh': {'tag': 'binding_protein', 'size': [[0, 71], [71, 100000]]},
                 'CAID23uh': {'tag': 'binding_protein', 'size': [[0, 71], [71, 100000]]},
                 'DBsh': {'tag': 'PDB', 'size': [[0, 71], [71, 100000]]}}
    _p = 'Data/'
    with open(f"Data/results/Figure_3/table.tsv", 'w') as fout:
        print(f"Data\tfrom\tto\tPercentage", file=fout)
        print(f"TR2008\t0\t30\t100%", file=fout)
        print(f"TR2017\t0\t30\t100%", file=fout)
        print(f"Data\tfrom\tto\tPercentage")
        for in_data in data_dict:
            tag = data_dict[in_data]['tag']
            af = aff_load3(f"{_p}af/{in_data}.af")
            aff_gen_counts(af)
            total = af['metadata']['counts']['tags_dict'][tag]['1']
            for rg in data_dict[in_data]['size']:
                af = aff_load3(f"{_p}af/{in_data}.af")
                aff_tag_size(af, tag, sz_range=rg)
                aff_gen_counts(af)
                cnt = af['metadata']['counts']['tags_dict'][tag]['1']
                print(f"{in_data}\t{rg[0]}\t{rg[1] - 1}\t{cnt / total:0.1%}", file=fout)
                print(f"{in_data}\t{rg[0]}\t{rg[1] - 1}\t{cnt / total:0.1%}")
