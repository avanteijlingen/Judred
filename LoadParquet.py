# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 03:29:47 2021

@author: avtei
"""

import pandas, sys, time
import numpy as np
import pyarrow.parquet as pq

Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

L = 6
fname = Num2Word[L].lower()+"peptides_normalized.parquet"

parquet_file = pq.ParquetFile(fname)

nRowGroups = parquet_file.num_row_groups

index_lower = 0
chunksize = -1
for i in range(nRowGroups):
    table = parquet_file.read_row_group(i)
    table = table.to_pandas()
    if chunksize == -1:
        chunksize = table.shape[0]
        index_upper = chunksize
    if chunksize > table.shape[0]: #sometimes there is a smaller one at the end
        index_upper = table.shape[0] + index_lower
    index = np.arange(index_lower, index_upper)
    table = table.set_index(index)
    print(table)
    time.sleep(1) # let the ram recover quickly
    index_lower += chunksize
    index_upper += chunksize

sys.exit()


Jparameters = pandas.read_parquet(fname)
if Jparameters.shape[0] != 20**L:
    print("Jparameters is missing rows", Jparameters.shape[0], "vs", 20**L)
    sys.exit()
print(Jparameters)
print("index type:", Jparameters.index.dtype)
if Jparameters.shape[0] == 0:
    sys.exit()
print(sys.getsizeof(Jparameters)/1024/1024, "MB")

Jparameters["SP2"] = Jparameters["SP2"].astype(np.float16)
Jparameters["NH2"] = Jparameters["NH2"].astype(np.float16)
Jparameters["S"] = Jparameters["S"].astype(np.float16)
Jparameters["Z"] = Jparameters["Z"].astype(np.float16)
Jparameters["OH"] = Jparameters["OH"].astype(np.float16)
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

    
a="""
Jparameters = Jparameters.astype(np.float16)
print(sys.getsizeof(Jparameters)/1024/1024, "MB")
"""
