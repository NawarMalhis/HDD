from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *


if __name__ == '__main__':
    t_head = 'Dataset\tTotalSeq\tTotalResidues\tSeq.WithSites\tNumberOfSites\tClass_0\tClass_1\tMasked_Residues'
    with open("Data/results/Tables/Table_S01_Data_Compositions.tsv", 'w') as fout:
        print(f"{t_head}", file=fout)
        for fl in fl_dict:
            af = aff_load3(in_file=f"{data_path}af/{fl}.af")
            aff_gen_counts(af)
            tag = fl_dict[fl]
            _c0 = af['metadata']['counts']['tags_dict'][tag]['0']
            _c1 = af['metadata']['counts']['tags_dict'][tag]['1']
            _cm = af['metadata']['counts']['tags_dict'][tag]['-']
            _seq_w_sites = af['metadata']['counts']['tags_dict'][tag]['seq']
            _sites = af['metadata']['counts']['tags_dict'][tag]['seg']
            total = _c0 + _cm + _c1
            print(f"{fl}\t{len(af['data'])}\t{total}\t{_seq_w_sites}\t{_sites}\t{_c0}\t{_c1}\t{_cm}", file=fout)

