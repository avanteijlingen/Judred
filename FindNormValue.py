# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 21:20:24 2021

@author: avtei
"""


import pandas, sys, os, time, math
import numpy as np




Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

L = 2



letters_1 = np.array(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
numbers = np.arange(0, len(letters_1), dtype=np.uint8)
features = ["SP2", "NH2", "MW", "S", "LogP WW", "Z", "MaxASA", "RotRatio", "Bulkiness", "OH"]
#sp2 carbons in side-chain only
SP2 =       np.array([0,    0,   1,   1,   6,   0,   3,   0,   0,   0,   0,   1,   0,   1,   1,   0,   0,   0,   8,   6], dtype=np.float32)
SP3 =       np.array([1,    1,   1,   2,   1,   0,   1,   4,   4,   4,   3,   1,   3,   2,   3,   1,   2,   3,   1,   1], dtype=np.float32)
NH2 =       np.array([0,    0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   0,   1,   2,   0,   0,   0,   0,   0], dtype=np.float32)
MW =        np.array([89.10, 121.16, 133.11, 147.13, 165.19, 75.07, 155.16, 131.18, 146.19, 131.18, 149.21, 132.12, 115.13, 146.15, 174.20, 105.09, 119.12, 117.15, 204.23, 181.19], dtype=np.float32)
S =         np.array([0,    1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0], dtype=np.float32)
charge =    np.array([0,    0,  -1,  -1,  0,    0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0], dtype=np.float32)
# ASP ARG and LYS are as charged side chains
Gwif = np.array([0.17, -0.24, 1.23, 2.02, -1.13, 0.01, 0.17, -0.31, 0.99, -0.56, -0.23, 0.42, 0.45, 0.58, 0.81, 0.13, 0.14, 0.07, -1.85, -0.94, ], dtype=np.float32) #kcal / mol
Gwoct = np.array([0.5, -0.02, 3.64, 3.63, -1.71, 1.15, 0.11, -1.12, 2.8, -1.25, -0.67, 0.85, 0.14, 0.77, 1.81, 0.46, 0.25, -0.46, -2.09, -0.71, ], dtype=np.float32) #kcal / mol
#Tien et al. 2013 (theory)
MaxASA =    np.array([129, 167, 193, 223, 240, 104, 224, 197, 236, 201, 224, 195, 159, 225, 274, 155, 172, 174, 285, 263], dtype=np.float32)
# Zimmerman J.M., Eliezer N., Simha R. J. Theor. Biol. 21:170-201(1968).
bulky =     np.array([11.50, 13.46, 11.68, 13.57, 19.80, 3.4, 13.69, 21.40, 15.71, 21.4, 16.25, 12.82, 17.43, 14.45, 14.28, 9.47, 15.77, 21.57, 21.67, 18.03], dtype=np.float32)
OH =        np.array([0,  0,   0,   0,    0,     0,    0,   0,   0,  0,   0,  0,  0,    0,   0,  1,  1,    0,   0,  1], dtype=np.float32)



SP2_min = ((min(SP2)*L)/2.0).astype(np.float32) 
SP2_max = ((max(SP2)*L)/2.0).astype(np.float32) 
polytryptophan_index = [18]*L
RotRatio_max = SP2[polytryptophan_index].sum() / SP3[polytryptophan_index].sum() 
RotRatio_max = np.float32(RotRatio_max/2.0)
NH2_max = ((max(NH2)*L)/2.0).astype(np.float32)
MW_min = (min(MW)*L).astype(np.float32)
MW_max = (max(MW)*L).astype(np.float32)
S_max = ((max(S)*L)/2.0).astype(np.float32)
Z_min = (min(charge)*L).astype(np.float32)
Z_max = (max(charge)*L).astype(np.float32)
polyasparticacid_index = [2]*L
LogP_WW_min = (Gwif[polyasparticacid_index] - Gwoct[polyasparticacid_index]).sum()
polyisoleucine_index = [7]*L
LogP_WW_max = (Gwif[polyisoleucine_index] - Gwoct[polyisoleucine_index]).sum()
MaxASA_min = (min(MaxASA)*L).astype(np.float32)
MaxASA_max = (max(MaxASA)*L).astype(np.float32)
bulky_min = (min(bulky)*L).astype(np.float32)
bulky_max = (max(bulky)*L).astype(np.float32)
OH_max = ((max(OH)*L)/2.0).astype(np.float32)



for LogP in np.arange(LogP_WW_min, LogP_WW_max, 0.2):
    print(LogP, end=" - ")
    LogP -= LogP_WW_min
    LogP /= ((LogP_WW_max - LogP_WW_min)/2.0)
    LogP -= 1.0
    LogP = LogP.astype(np.float32)
    print(LogP)

print("Z")
for Z in [-1, 0, 1]: #np.arange(Z_min, Z_max, 0.2):
    print(Z, end=" - ")
    Z -= Z_min
    Z /= ((Z_max - Z_min)/2.0)
    Z -= 1.0
    Z = Z.astype(np.float32)
    print(Z)

