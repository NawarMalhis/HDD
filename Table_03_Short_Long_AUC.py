from param import *
import sys
if aff_path not in sys.path:
    sys.path.append(aff_path)
    
from annotated_fasta import *
from annotated_fasta_metrics import *
from annotated_fasta_CAID import *


if __name__ == '__main__':
    prd_list = ['AlphaFold-binding', 'ANCHOR-2', 'CNN_C1u', 'CNN_C23u', 'CNN_DBs', 'CNN_TR08u', 'DeepDISObind-protein',
                'DeepDRPBind-protein', 'DisoRDPbind-protein', 'DRPBind-protein', 'fMoRFpred', 'OPAL', 'MoRFchibi',
                'MoRFchibi-light', 'MoRFchibi-web']
    files_dict = {'CAID1uh': 'binding_protein', 'CAID23uh': 'binding_protein', 'DBsh': 'PDB'}

    for sl in ['_short', '_long']:
        auc_dict = {}
        for fl in files_dict:
            prd_used = dataset_prd_dict[fl]
            tag = files_dict[fl]
            af = aff_load3(in_file=f"{data_path}af/{fl}.af")
            aff_load_caid_scores(af, scores_path=f"{data_path}scores/", prd_list=prd_used, merged=False,
                                 remove_missing_scores=False)
            if sl == '_long':
                aff_tag_size(af, tag, sz_range=[71, 10000])  # mask out 'tag's <71 AAs
            else:
                aff_tag_size(af, tag, sz_range=[0, 71])  # mask out 'tag's >=71 AAs
            auc_dict[fl], _ = aff_roc(af, tag=tag, prd_list=prd_used, display=False)

        with open(f'Data/results/Tables/Table_3{sl}.tsv', 'w') as fout:
            print("Predictor\tGroup\tCAID1uh\tCAID2&3uh\tDBsh", file=fout)
            for prd in prd_list:
                print(f"{prd}", end=':\t', file=fout)
                for fl in files_dict:
                    if prd in auc_dict[fl]:
                        print(f"|{auc_dict[fl][prd]:.3f}", end='\t', file=fout)
                    else:
                        print(end='|\t', file=fout)
                print(file=fout)

