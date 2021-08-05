# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 03:29:47 2021

@author: avtei
"""

import pandas, sys
import numpy as np

Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

L = 3

Jparameters = pandas.read_hdf(Num2Word[L].lower()+"peptides_normalized.hdf5")

print(Jparameters)

if Jparameters.shape[0] != 20**L:
    print("Jparameters is missing rows", Jparameters.shape[0], "vs", 20**L)
    sys.exit()
#print(Jparameters)
print("index type:", Jparameters.index.dtype)
if Jparameters.shape[0] == 0:
    sys.exit()
print(sys.getsizeof(Jparameters)/1024/1024, "MB")

# =============================================================================
# Jparameters["SP2"] = Jparameters["SP2"].astype(np.float16)
# Jparameters["NH2"] = Jparameters["NH2"].astype(np.float16)
# Jparameters["S"] = Jparameters["S"].astype(np.float16)
# Jparameters["Z"] = Jparameters["Z"].astype(np.float16)
# Jparameters["OH"] = Jparameters["OH"].astype(np.float16)
# =============================================================================
print(sys.getsizeof(Jparameters)/1024/1024, "MB")

for col in Jparameters.columns:
    print(col, Jparameters[col].max(), Jparameters[col].min(), Jparameters[col].dtype)

print(np.float32)
for col in ["MW", "S"]:
    print(col, "unique values:",  np.unique(Jparameters[col]).shape[0])
Jparameters = Jparameters.astype(np.float16)
print(np.float16)
for col in ["MW", "S"]:
    print(col, "unique values:",  np.unique(Jparameters[col]).shape[0])
