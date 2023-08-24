# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 13:08:32 2023

@author: Alex
"""

import pandas, sys, os, time, math, re
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from datetime import datetime
import orca_parser

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


try:
    import cupy as cp
    mempool = cp.get_default_memory_pool()
    mempool.set_limit(size=3*1024**3) # 3 GB
    use_gpu = True
except:
    use_gpu = False


Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

# =============================================================================
# L = int(sys.argv[1])
# steps = np.array([20**i for i in range(L-1, -1, -1)], dtype=np.uint64)
# 
# fname = Num2Word[L].lower()+"peptides_normalized.parquet"
# #fname_memmap = Num2Word[L].lower()+"peptides_indexes.memmap"
# =============================================================================
# Nfex = 'Nfer'/'Nfes',
letters_1 = np.array(['Na', 'Nab', 'Nd', 'Ne', 'Nf', 'Nfe', 'NfeBr', 'NfeCl', "Nfex", 'Nfn', 'Nfnap', 'Nfp', 'Ni', 'Nk', 'Nke', 'Nkeqm', 'Nl', 'Nm', 'NmO', 'Nn', 'Nq', 'Nr', 'Ns', 'Nse', 'Nt', 'Nv', 'Nw', 'Nwe', 'Ny'])

for peptoid in letters_1:
    if peptoid == "Nfex":
        peptoid_folder = "Nfes"
    else:
        peptoid_folder = peptoid
    
    for state in ['gas', 'octanol', 'water']:
        folder = f"peptoid_data/{peptoid_folder}/{state}"
        for file in os.listdir(folder):
            path = f"{folder}/{file}"
            op = orca_parser.ORCAParse(path)
            print(op.valid)
    
    
    break