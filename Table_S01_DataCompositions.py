#!/usr/bin/env python3
"""
Supplementary Table S01 — Dataset composition statistics
(Total sequences, residues, binding sites, class counts, etc.)

Author: Nawar Malhis
Refined using Grok
The University of British Columbia, 2026
"""

from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *


if __name__ == '__main__':
    fl_dict = {
        'DisProt': 'binding_protein',
        'CAID1u': 'binding_protein',
        'CAID1uh': 'binding_protein',
        'CAID23u': 'binding_protein',
        'CAID23uh': 'binding_protein',
        'DBs': 'PDB',
        'DBsh': 'PDB',
        'TR2008u': 'PDB',
        'VA': 'binding_protein',
        'VA_DisProt': 'binding_protein',
        'VA_PDB': 'binding_protein',
    }

    t_head = (
        "Dataset\tTotalSeq\tSeq.WithSites\tTotalResidues\t"
        "NumberOfSites\tClass_0\tClass_1\tMasked_Residues"
    )

    out_file = "Data/results/Tables/Table_S01_Data_Compositions.tsv"

    with open(out_file, 'w') as fout:
        print(t_head, file=fout)

        for fl in fl_dict:
            af = aff_load3(in_file=f"{data_path}af/{fl}.af")
            aff_gen_counts(af)

            tag = fl_dict[fl]

            c0 = af['metadata']['counts']['tags_dict'][tag]['0']
            c1 = af['metadata']['counts']['tags_dict'][tag]['1']
            cm = af['metadata']['counts']['tags_dict'][tag]['-']
            seq_w_sites = af['metadata']['counts']['tags_dict'][tag]['seq']
            num_sites = af['metadata']['counts']['tags_dict'][tag]['seg']
            total_res = c0 + cm + c1

            print(
                f"{fl}\t"
                f"{len(af['data'])}\t"
                f"{seq_w_sites}\t"
                f"{total_res}\t"
                f"{num_sites}\t"
                f"{c0}\t"
                f"{c1}\t"
                f"{cm}",
                file=fout
            )

    print(f"Table saved to: {out_file}")