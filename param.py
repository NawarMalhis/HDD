
aff_path = '/home/nmalhis/tmp/AFF/'
data_path = "Data/"

fl_dict = {'CAID1u': 'binding_protein', 'CAID1uh': 'binding_protein', 'CAID23u': 'binding_protein',
           'CAID23uh': 'binding_protein', 'DBs': 'PDB', 'DBsh': 'PDB', 'TR2008u': 'PDB', 'VA': 'binding_protein',
           'VA_DisProt': 'binding_protein', 'VA_PDB': 'binding_protein'}

dataset_prd_dict = {'DBsh': ['AlphaFold-binding', 'ANCHOR-2', 'CNN_C1u', 'CNN_C23u', 'CNN_TR08u',
                             'DeepDISObind-protein', 'DeepDRPBind-protein', 'DisoRDPbind-protein',
                             'DRPBind-protein', 'fMoRFpred', 'MoRFchibi', 'MoRFchibi-light', 'MoRFchibi-web',
                             'OPAL'],
                    'CAID23uh': ['AlphaFold-binding', 'ANCHOR-2', 'CNN_C1u', 'CNN_DBs', 'CNN_TR08u',
                                 'DeepDISObind-protein',
                                 'DeepDRPBind-protein', 'DisoRDPbind-protein', 'DRPBind-protein', 'fMoRFpred',
                                 'MoRFchibi',
                                 'MoRFchibi-light', 'MoRFchibi-web', 'OPAL'],
                    'CAID1uh': ['ANCHOR-2', 'CNN_C23u', 'CNN_DBs', 'CNN_TR08u', 'DisoRDPbind-protein', 'fMoRFpred',
                                'MoRFchibi', 'MoRFchibi-light', 'MoRFchibi-web', 'OPAL']}

prd_dict = {'ANCHOR-2': {'weight': 2, 'style': '-', 'color': 'lime'},
           'MoRFchibi': {'weight': 2, 'style': ':', 'color': 'g'},
           'MoRFchibi-light': {'weight': 2, 'style': '--', 'color': 'g'},
           'MoRFchibi-web': {'weight': 2, 'style': '-', 'color': 'g'},
           'OPAL': {'weight': 2, 'style': '-', 'color': 'b'},
           'DisoRDPbind-protein': {'weight': 2, 'style': '-', 'color': 'brown'},
           'fMoRFpred': {'weight': 2, 'style': ':', 'color': 'brown'},
           'DeepDISObind-protein': {'weight': 2, 'style': '-', 'color': 'blueviolet'},
           'AlphaFold-binding': {'weight': 2, 'style': '-', 'color': 'coral'},
           'DeepDRPBind-protein': {'weight': 2, 'style': '-', 'color': 'tomato'},
           'DRPBind-protein': {'weight': 2, 'style': '-', 'color': 'y'},
            'CNN_C1u': {'weight': 2, 'style': '-', 'color': 'r'},
            'CNN_C23u': {'weight': 2, 'style': '--', 'color': 'r'},
            'CNN_DBs': {'weight': 2, 'style': ':', 'color': 'r'},
            'CNN_TR08u': {'weight': 2, 'style': '-.', 'color': 'r'},
            }

