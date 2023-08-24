# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 13:08:32 2023

@author: Alex
"""

import pandas, sys, os, time, math, re, glob
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

peptoids = pandas.DataFrame(index=['Na', 'Nab', 'Nd', 'Ne', 'Nf', 'Nfe', 'NfeBr', 'NfeCl', "Nfex", 'Nfn', 'Nfnap', 'Nfp', 'Ni', 'Nk', 'Nke', 'Nkeqm', 'Nl', 'Nm', 'NmO', 'Nn', 'Nq', 'Nr', 'Ns', 'Nse', 'Nt', 'Nv', 'Nw', 'Nwe', 'Ny'],
                            columns=["Gwif", "Gwoct"])
# =============================================================================
# L = int(sys.argv[1])
# steps = np.array([20**i for i in range(L-1, -1, -1)], dtype=np.uint64)
# 
# fname = Num2Word[L].lower()+"peptides_normalized.parquet"
# #fname_memmap = Num2Word[L].lower()+"peptides_indexes.memmap"
# =============================================================================
# Nfex = 'Nfer'/'Nfes',

for peptoid in peptoids.index:
    if peptoid == "Nfex":
        peptoid_folder = "Nfes"
    else:
        peptoid_folder = peptoid
    
    for i in range(3):
        gas = orca_parser.ORCAParse(glob.glob(f"{peptoid_folder}/gas/*.out")[i])
        if gas.valid:
            break
    for i in range(3):
        octanol = orca_parser.ORCAParse(glob.glob(f"{peptoid_folder}/octanol/*.out")[i])
        if octanol.valid:
            break
    for i in range(3):
        water = orca_parser.ORCAParse(glob.glob(f"{peptoid_folder}/water/*.out")[i])
        if water.valid:
            break

    for o in [gas, octanol, water]:
        print(o.fname, o.valid)
        o.parse_free_energy()
    

    peptoids.at[peptoid, "Gwif"] = water.Gibbs - gas.Gibbs
    peptoids.at[peptoid, "Gwoct"] = octanol.Gibbs - gas.Gibbs
    
# =============================================================================
#     for state in ['gas', 'octanol', 'water']:
#         folder = f"{peptoid_folder}/{state}"
#         for file in os.listdir(folder):
#             path = f"{folder}/{file}"
#             op = orca_parser.ORCAParse(path)
#             print(op.valid)
#             if not op.valid:
#                 os.remove(path)
#     break
# =============================================================================
peptoids["Solubility Scale"] = peptoids["Gwif"] - peptoids["Gwoct"]

ASA = pandas.read_csv("ASA.csv", index_col=0)
for index in ASA.index:
    peptoids.at[index, "ASA"] = ASA.at[index, "ASA"]


x = pandas.read_csv("judtoid.csv")
x = x[x["State"] == "water"]
x.index = x["Species"]
x = x.sort_index()
peptoids = peptoids.sort_index() 
for col in ["S", "OH", "Charge", "Mass (AMU)", "NH2", "SP2", "SP3"]:
    peptoids[col] = x[col]
    peptoids.at["Nfex", col] = x.at["Nfes", col]
    

print(peptoids.sort_values("Solubility Scale"))


peptoids.to_csv("Peptoid_Data.csv")






