# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 23:36:37 2021

@author: avtei
"""

a="""
import sqlite3 as sl

con = sl.connect("tripeptides.db")


con.close()
"""
import pandas, itertools, sys, os, h5py, time
import numpy as np
from dask.distributed import Client, progress
import dask.dataframe as dd
import dask.array as da

st = time.time()

try:
    progress() #Raises value error if client isn't already running
except ValueError:
    client = Client(processes=False, threads_per_worker=1, n_workers=2, memory_limit='1GB') 



Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

features = ["SP2", "NH2", "MW", "S", "LogP WW", "Z", "MaxASA", "RotRatio", "Bulkiness", "OH"]
#sp2 carbons in side-chain only
SP2 =       np.array([0,    0,   1,   1,   6,   0,   3,   0,   0,   0,   0,   1,   0,   1,   1,   0,   0,   0,   8,   6])
SP3 =       np.array([1,    1,   1,   2,   1,   0,   1,   4,   4,   4,   3,   1,   3,   2,   3,   1,   2,   3,   1,   1])
NH2 =       np.array([0,    0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   0,   1,   2,   0,   0,   0,   0,   0])
MW =        np.array([89.10, 121.16, 133.11, 147.13, 165.19, 75.07, 155.16, 131.18, 146.19, 131.18, 149.21, 132.12, 115.13, 146.15, 174.20, 105.09, 119.12, 117.15, 204.23, 181.19])
S =         np.array([0,    1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0])
charge =    np.array([0,    0,  -1,  -1,  0,    0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0])
# ASP ARG and LYS are as charged side chains
Gwif = np.array([0.17, -0.24, 1.23, 2.02, -1.13, 0.01, 0.17, -0.31, 0.99, -0.56, -0.23, 0.42, 0.45, 0.58, 0.81, 0.13, 0.14, 0.07, -1.85, -0.94, ]) #kcal / mol
Gwoct = np.array([0.5, -0.02, 3.64, 3.63, -1.71, 1.15, 0.11, -1.12, 2.8, -1.25, -0.67, 0.85, 0.14, 0.77, 1.81, 0.46, 0.25, -0.46, -2.09, -0.71, ]) #kcal / mol
#Tien et al. 2013 (theory)
MaxASA =    np.array([129, 167, 193, 223, 240, 104, 224, 197, 236, 201, 224, 195, 159, 225, 274, 155, 172, 174, 285, 263])
# Zimmerman J.M., Eliezer N., Simha R. J. Theor. Biol. 21:170-201(1968).
bulky =     np.array([11.50, 13.46, 11.68, 13.57, 19.80, 3.4, 13.69, 21.40, 15.71, 21.4, 16.25, 12.82, 17.43, 14.45, 14.28, 9.47, 15.77, 21.57, 21.67, 18.03])
OH =        np.array([0,  0,   0,   0,    0,     0,    0,   0,   0,  0,   0,  0,  0,    0,   0,  1,  1,    0,   0,  1])
#nBase =     np.array([0,   0,  0,   0,    0,     0,    0,   0,   1,  0,   0,  0,  0,    0,   3,    0,  0,  0,   0,  0])


L = 3

letters_1 = np.array(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
numbers = np.arange(0, len(letters_1))


# Since it takes quite a whole to generate the peptide sequences we'll save them incase it crashes or stops for any reason it can be restarted
if os.path.exists(Num2Word[L]+"peptide_names.npy"):
    peptides = np.load(Num2Word[L]+"peptide_names.npy")
    peptide_numbers = np.load(Num2Word[L]+"peptide_numbers.npy")
else:
    nt = time.time()
    peptides = np.array(list(itertools.product(letters_1, repeat = L)))
    peptides = da.array(["".join(LETTERS) for LETTERS in peptides])
    peptides = peptides.compute()
    peptide_numbers = np.array(list(itertools.product(numbers, repeat = L)))
    np.save(Num2Word[L]+"peptide_numbers.npy", peptide_numbers)
    np.save(Num2Word[L]+"peptide_names.npy", peptides)
    print(Num2Word[L]+"peptides names generated in", round(time.time()-nt, 3), "s")


SP2_data = da.array(SP2[peptide_numbers]).sum(axis=1).compute()
NH2_data = da.array(NH2[peptide_numbers]).sum(axis=1).compute()
MW_data = da.array(MW[peptide_numbers]).sum(axis=1).compute()
S_data = da.array(S[peptide_numbers]).sum(axis=1).compute()

LogPWW_data = da.array(Gwif[peptide_numbers]) - da.array(Gwoct[peptide_numbers]).compute()
LogPWW_data = LogPWW_data.sum(axis=1).compute()

Z_data = da.array(charge[peptide_numbers]).sum(axis=1).compute()
MaxASA_data = da.array(MaxASA[peptide_numbers]).sum(axis=1).compute()

RotRatio_data = da.array(SP3[peptide_numbers]).sum(axis=1)
RotRatio_data = da.array(SP2_data / RotRatio_data).compute()

bulky_data = da.array(bulky[peptide_numbers]).sum(axis=1).compute()
OH_data = da.array(OH[peptide_numbers]).sum(axis=1).compute()

data = da.vstack((SP2_data, NH2_data, MW_data, S_data, LogPWW_data, Z_data,
                  MaxASA_data, RotRatio_data, bulky_data, OH_data)).T
#npdata = data.compute()
#print(npdata)


Jparameters = dd.from_dask_array(data, columns=features).compute()
Jparameters = Jparameters.set_index(peptides)
Jparameters.to_parquet(Num2Word[L]+"peptides.parquet")


print(Num2Word[L]+"peptides dataset done in", round(time.time()-st, 3), "s")

for file in [Num2Word[L]+"peptide_names.npy", Num2Word[L]+"peptide_numbers.npy"]:
    if os.path.exists(file):
        os.remove(file)