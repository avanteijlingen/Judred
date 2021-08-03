# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 21:20:24 2021

@author: avtei
"""


import pandas, sys, os, time
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

def index2pep(index, Length):
    size = 20**Length
    if index >= size:
        print("No peptide this large in the dataset")
        return None
    letters_1 = list("ACDEFGHIKLMNPQRSTVWY")
    step = int(size/20)
    above_size = int(size/20)
    steps = np.array([i*step for i in range(20)])
    solution = []
    letter_i = np.where(steps <= index)[0][-1]
    letter = letters_1[letter_i]
    solution.append(letter)
    while len(solution) < Length:
        index = index%step
        step = int(step/20)
        steps = np.array([i*step for i in range(20)])
        if index <= above_size/20 and len(solution) < Length-1:
            letter = "A"
            #print(index <= above_size/20)
        else:
            letter_i = np.where(steps <= index)[0][-1]
            letter = letters_1[letter_i]
            #print(letter_i, letter)
        above_size = int(above_size/20)
        solution.append(letter)
    return "".join(solution)

letters_1 = np.array(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
numbers = np.arange(0, len(letters_1), dtype=np.uint8)
features = ["SP2", "NH2", "MW", "S", "LogP WW", "Z", "MaxASA", "RotRatio", "Bulkiness", "OH"]
#sp2 carbons in side-chain only
SP2 =       np.array([0,    0,   1,   1,   6,   0,   3,   0,   0,   0,   0,   1,   0,   1,   1,   0,   0,   0,   8,   6], dtype=np.int8)
SP3 =       np.array([1,    1,   1,   2,   1,   0,   1,   4,   4,   4,   3,   1,   3,   2,   3,   1,   2,   3,   1,   1], dtype=np.int8)
NH2 =       np.array([0,    0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   0,   1,   2,   0,   0,   0,   0,   0], dtype=np.int8)
MW =        np.array([89.10, 121.16, 133.11, 147.13, 165.19, 75.07, 155.16, 131.18, 146.19, 131.18, 149.21, 132.12, 115.13, 146.15, 174.20, 105.09, 119.12, 117.15, 204.23, 181.19], dtype=np.float16)
S =         np.array([0,    1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0], dtype=np.int8)
charge =    np.array([0,    0,  -1,  -1,  0,    0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0], dtype=np.int8)
# ASP ARG and LYS are as charged side chains
Gwif = np.array([0.17, -0.24, 1.23, 2.02, -1.13, 0.01, 0.17, -0.31, 0.99, -0.56, -0.23, 0.42, 0.45, 0.58, 0.81, 0.13, 0.14, 0.07, -1.85, -0.94, ], dtype=np.float16) #kcal / mol
Gwoct = np.array([0.5, -0.02, 3.64, 3.63, -1.71, 1.15, 0.11, -1.12, 2.8, -1.25, -0.67, 0.85, 0.14, 0.77, 1.81, 0.46, 0.25, -0.46, -2.09, -0.71, ], dtype=np.float16) #kcal / mol
#Tien et al. 2013 (theory)
MaxASA =    np.array([129, 167, 193, 223, 240, 104, 224, 197, 236, 201, 224, 195, 159, 225, 274, 155, 172, 174, 285, 263], dtype=np.int16)
# Zimmerman J.M., Eliezer N., Simha R. J. Theor. Biol. 21:170-201(1968).
bulky =     np.array([11.50, 13.46, 11.68, 13.57, 19.80, 3.4, 13.69, 21.40, 15.71, 21.4, 16.25, 12.82, 17.43, 14.45, 14.28, 9.47, 15.77, 21.57, 21.67, 18.03], dtype=np.float16)
OH =        np.array([0,  0,   0,   0,    0,     0,    0,   0,   0,  0,   0,  0,  0,    0,   0,  1,  1,    0,   0,  1], dtype=np.int8)
    
L = int(sys.argv[1])
a = np.indices((len(numbers),) * L, dtype=np.int8)
b = np.rollaxis(a, 0, L + 1)
c = b.reshape(-1, L)
print("c:", c.nbytes/1024/1024, "MB")


fname = Num2Word[L].lower()+"peptides.parquet"
a="""
index = [(20**L)-2, (20**L)-1]
peptide_numbers = c[index]
print(index, peptide_numbers)
print(index2pep(index[0], L), index2pep(index[1], L))
print("OH:", OH[peptide_numbers].sum(axis=1).astype(np.uint8))
print("OH:", OH[peptide_numbers].sum(axis=1).astype(np.uint8))
#"""
#a="""
chunksize = 1000000
y = np.zeros((chunksize, 10), dtype=np.float32)
pd_table = pandas.DataFrame(y, columns=features, index=np.arange(0,chunksize))
pd_table["SP2"] = pd_table["SP2"].astype(np.uint8)
pd_table["NH2"] = pd_table["NH2"].astype(np.uint8)
pd_table["MW"] = pd_table["MW"].astype(np.float32)
pd_table["S"] = pd_table["S"].astype(np.uint8)
pd_table["LogP WW"] = pd_table["LogP WW"].astype(np.float32)
pd_table["Z"] = pd_table["Z"].astype(np.int8)
pd_table["MaxASA"] = pd_table["MaxASA"].astype(np.uint16)
pd_table["RotRatio"] = pd_table["RotRatio"].astype(np.float32)
pd_table["Bulkiness"] = pd_table["Bulkiness"].astype(np.float32)
pd_table["OH"] = pd_table["OH"].astype(np.uint8)
table = pa.Table.from_pandas(pd_table, preserve_index=False)

with pq.ParquetWriter(fname, table.schema) as writer:
    for i in range(int((20**L)/chunksize)):
        pd_table = pandas.DataFrame(y, columns=features, index=np.arange(0,chunksize))
        
        peptide_numbers = c[i*chunksize:(i+1)*chunksize]
        pd_table["SP2"] = SP2[peptide_numbers].sum(axis=1).astype(np.uint8)
        pd_table["NH2"] = NH2[peptide_numbers].sum(axis=1).astype(np.uint8)
        pd_table["MW"] = MW[peptide_numbers].sum(axis=1).astype(np.float32)
        pd_table["S"] = S[peptide_numbers].sum(axis=1).astype(np.uint8)
        pd_table["LogP WW"] = (Gwif[peptide_numbers] - Gwoct[peptide_numbers]).sum(axis=1).astype(np.float32)
        pd_table["Z"] = charge[peptide_numbers].sum(axis=1).astype(np.int8)
        pd_table["MaxASA"] = MaxASA[peptide_numbers].sum(axis=1).astype(np.uint16)
        pd_table["RotRatio"] = (pd_table["SP2"]/(SP3[peptide_numbers].sum(axis=1))).astype(np.float32)
        pd_table["Bulkiness"] = bulky[peptide_numbers].sum(axis=1).astype(np.float32)
        pd_table["OH"] = OH[peptide_numbers].sum(axis=1).astype(np.uint8)
        
        table = pa.Table.from_pandas(pd_table, preserve_index=False)
        writer.write_table(table)
        
        print("{:.3e}".format(chunksize*(i+1)), "/", "{:.3e}".format(20**L))

#"""