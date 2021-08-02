# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 23:36:37 2021

@author: avtei
"""

import gc
import pandas, itertools, sys, os, h5py, time, gc
import numpy as np
from dask.distributed import Client, progress
import dask.dataframe as dd
import dask.array as da
import pyarrow as pa

st = time.time()

try:
    progress() #Raises value error if client isn't already running
except ValueError:
    #client = Client('127.0.0.1:8787')
    client = Client(processes=False, threads_per_worker=1, n_workers=1, memory_limit='1000GB', silence_logs='error') 
    print(client)



Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

data_type = np.float16 #Parquet doesn't currently support float16, so we have to save as float32 then convert when we load it in again

features = ["SP2", "NH2", "MW", "S", "LogP WW", "Z", "MaxASA", "RotRatio", "Bulkiness", "OH"]
#sp2 carbons in side-chain only
SP2 =       np.array([0,    0,   1,   1,   6,   0,   3,   0,   0,   0,   0,   1,   0,   1,   1,   0,   0,   0,   8,   6], dtype=data_type)
SP3 =       np.array([1,    1,   1,   2,   1,   0,   1,   4,   4,   4,   3,   1,   3,   2,   3,   1,   2,   3,   1,   1], dtype=data_type)
NH2 =       np.array([0,    0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   0,   1,   2,   0,   0,   0,   0,   0], dtype=data_type)
MW =        np.array([89.10, 121.16, 133.11, 147.13, 165.19, 75.07, 155.16, 131.18, 146.19, 131.18, 149.21, 132.12, 115.13, 146.15, 174.20, 105.09, 119.12, 117.15, 204.23, 181.19], dtype=data_type)
S =         np.array([0,    1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0], dtype=data_type)
charge =    np.array([0,    0,  -1,  -1,  0,    0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0], dtype=data_type)
# ASP ARG and LYS are as charged side chains
Gwif = np.array([0.17, -0.24, 1.23, 2.02, -1.13, 0.01, 0.17, -0.31, 0.99, -0.56, -0.23, 0.42, 0.45, 0.58, 0.81, 0.13, 0.14, 0.07, -1.85, -0.94, ], dtype=data_type) #kcal / mol
Gwoct = np.array([0.5, -0.02, 3.64, 3.63, -1.71, 1.15, 0.11, -1.12, 2.8, -1.25, -0.67, 0.85, 0.14, 0.77, 1.81, 0.46, 0.25, -0.46, -2.09, -0.71, ], dtype=data_type) #kcal / mol
#Tien et al. 2013 (theory)
MaxASA =    np.array([129, 167, 193, 223, 240, 104, 224, 197, 236, 201, 224, 195, 159, 225, 274, 155, 172, 174, 285, 263], dtype=data_type)
# Zimmerman J.M., Eliezer N., Simha R. J. Theor. Biol. 21:170-201(1968).
bulky =     np.array([11.50, 13.46, 11.68, 13.57, 19.80, 3.4, 13.69, 21.40, 15.71, 21.4, 16.25, 12.82, 17.43, 14.45, 14.28, 9.47, 15.77, 21.57, 21.67, 18.03], dtype=data_type)
OH =        np.array([0,  0,   0,   0,    0,     0,    0,   0,   0,  0,   0,  0,  0,    0,   0,  1,  1,    0,   0,  1], dtype=data_type)
#nBase =     np.array([0,   0,  0,   0,    0,     0,    0,   0,   1,  0,   0,  0,  0,    0,   3,    0,  0,  0,   0,  0])


L = int(sys.argv[1])

letters_1 = np.array(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
numbers = da.arange(0, len(letters_1), dtype=np.int8)


# Since it takes quite a whole to generate the peptide sequences we'll save them incase it crashes or stops for any reason it can be restarted
if os.path.exists(Num2Word[L]+"peptide_names.npy") and 1 ==0:
    peptides = np.load(Num2Word[L]+"peptide_names.npy")
    peptide_numbers = np.load(Num2Word[L]+"peptide_numbers.npy")
else:
    nt = time.time()
    #peptides = da.array(["".join(LETTERS) for LETTERS in list(itertools.product(letters_1, repeat = L))], dtype=dd.Index)
    #peptide_numbers = da.array(list(itertools.product(numbers, repeat = L)))
    
    #np.array(list(itertools.product(some_list, repeat=some_length)))
    
    a = da.indices((len(numbers),) * L, dtype=np.int8)
    b = da.rollaxis(a, 0, L + 1)
    c = b.reshape(-1, L)

    print(Num2Word[L]+"peptides names/numbers generated in", round(time.time()-nt, 3), "s")
#print("Peptide numbers in RAM:", x.nbytes/1024/1024, "MB")

#data = da.zeros((20**L,10), dtype=np.float16)
#dask_df = dd.from_dask_array(data)

#dask_df[0] = [0]*(20**L)
    
# Do biggest tasks first so the final 1-D array is collected first to keep more memory free later down the line
LogPWW_data = da.array(da.array(Gwif[c]) - da.array(Gwoct[c]))
LogPWW_data = LogPWW_data.sum(axis=1)
LogPWW_data = da.array(LogPWW_data.compute())
print("LogP_WW generated")

SP2_data = da.array(SP2[c]).sum(axis=1)
SP2_data = da.array(SP2_data.compute()) # by keeping it inside a da.array memory will not overflow
print("LogP_WW generated")

RotRatio_data = da.array(SP3[c]).sum(axis=1)
RotRatio_data = da.array(SP2_data / RotRatio_data)
RotRatio_data = RotRatio_data.compute()
np.nan_to_num(RotRatio_data, copy=False, nan=0.0) # fillna is for polyglycine which results in nan cus its devide by 0
RotRatio_data = da.array(RotRatio_data) 

print("RotRatio generated")

NH2_data = da.array(NH2[c]).sum(axis=1)
NH2_data = da.array(NH2_data.compute())
print("NH2 generated")

MW_data = da.array(MW[c]).sum(axis=1)
MW_data = da.array(MW_data.compute())
print("MW generated")

S_data = da.array(S[c]).sum(axis=1)
S_data = da.array(S_data.compute())
print("S generated")

Z_data = da.array(charge[c]).sum(axis=1)
Z_data = da.array(Z_data.compute())
print("Z generated")

MaxASA_data = da.array(MaxASA[c]).sum(axis=1)
MaxASA_data = da.array(MaxASA_data.compute())
print("MaxASA generated")

bulky_data = da.array(bulky[c]).sum(axis=1)
bulky_data = da.array(bulky_data.compute())
print("bulky generated")

OH_data = da.array(OH[c]).sum(axis=1)
OH_data = da.array(OH_data.compute())
print("OH generated")




data = da.vstack((SP2_data, NH2_data, MW_data, S_data, LogPWW_data, Z_data,
                  MaxASA_data, RotRatio_data, bulky_data, OH_data)).T

Jparameters = dd.from_dask_array(data, columns=features).compute()
del LogPWW_data
del SP2_data
del RotRatio_data
del NH2_data
del MW_data
del S_data
del Z_data
del MaxASA_data
del bulky_data
del OH_data

gc.collect()

Jparameters = Jparameters.astype(np.float32)
peptides = da.array(["".join(x) for x in letters_1[c]]).compute()
Jparameters = Jparameters.set_index(peptides)
Jparameters.to_parquet(Num2Word[L].lower()+"peptides.parquet")

#Jparameters = dd.from_dask_array(data, columns=features).compute()
#Jparameters = Jparameters.set_index(peptides)
#Jparameters.to_parquet(Num2Word[L]+"peptides.parquet")



print(Num2Word[L]+"peptides dataset done in", round(time.time()-st, 3), "s")


a="""
for file in [Num2Word[L]+"peptide_names.npy", Num2Word[L]+"peptide_numbers.npy"]:
    if os.path.exists(file):
        os.remove(file)
"""
