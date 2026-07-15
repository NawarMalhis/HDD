from annotated_fasta_CAID import aff_load_caid_scores
from annotated_fasta import aff_load3, aff_save3, aff_tag_size
from annotated_fasta_metrics import aff_roc, aff_violin_plot
import numpy as np
from param import *


if __name__ == '__main__':
    _cnn = False
    _p = 'Data/'
    for in_data in ['DBsh', 'CAID1uh', 'CAID23uh']:
        tag = 'binding_protein'
        if 'DBs' in in_data:
            tag = 'PDB'
        for sl in ['']:  # sl stands for '_short', '_long']:
            target_data = f"{in_data}{sl}"
            tools_list = ['ANCHOR-2', 'MoRFchibi', 'OPAL', 'MoRFchibi-light', 'MoRFchibi-web', 'DisoRDPbind-protein',
                          'fMoRFpred']

            af = aff_load3(f"{_p}af/{target_data}.af")  # _UniParc  _extra
            aff_tag_size(af, tag, sz_range=[71, 10000])
            aff_save3(af, f"{_p}af/{target_data}_l71.af")




# aff_tag_size(af, tag, sz_range