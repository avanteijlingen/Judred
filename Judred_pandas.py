# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 00:15:57 2021

@author: avtei
"""


import numpy as np
import sys, time, pandas
import peptideutils as pu


print("THIS NEW/SPECIAL VERSION INCLUDES HELICAL PENALTY")

DoMinMaxScaling = False
pd_table = pandas.DataFrame()
# CPU

letters_1 = np.array(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
numbers = np.arange(0, len(letters_1), dtype=np.uint8)
features = ["SP2", "NH2", "MW", "S", "LogP WW", "Z", "MaxASA", "RotRatio", "Bulkiness", "OH", "pI"]
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
pI =        np.array([6.11,5.15,2.98,3.08,5.76,6.06, 7.64, 6.04, 9.47, 6.04, 5.71, 5.43, 6.30, 5.65, 11.5,5.07, 5.60,6.02,5.88,5.63], dtype=np.float32)

#Helical_penalty
#https://doi.org/10.1016%2Fs0006-3495%2898%2977529-0
Helical_penalty = np.array([0.00,0.68,0.69,0.40,0.54,1.0,0.61,0.41,0.26,0.21,0.24,0.65,3.16,0.39,0.21,0.50,0.66,0.61,0.49,0.53]) # kcal/mol

L = 2

if DoMinMaxScaling:
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
    pI_min = (min(pI)).astype(np.float32)
    pI_max = (max(pI)).astype(np.float32)
    Helical_penalty_max = (max(Helical_penalty)).astype(np.float32)
    Helical_penalty_min = (min(Helical_penalty)).astype(np.float32)


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
        #print("index:", index)
        step = int(step/20)
        steps = np.array([(i*step) for i in range(20)])
        #print(steps)
        if index < above_size/20 and len(solution) < Length-1:
            letter = "A"
            #print("index < above_size/20")
        else:
            letter_i = np.where(steps <= index)[0][-1]
            letter = letters_1[letter_i]
            #print(letter_i, letter)
        above_size = int(above_size/20)
        solution.append(letter)
    return "".join(solution)




peptides = pu.GenerateDatasetIndex(2)
peptide_numbers = np.zeros((len(peptides), L), dtype=np.uint8)
for i,pep in enumerate(peptides):
    for j,letter in enumerate(pep):
        peptide_numbers[i][j] = pu.pep2index(letter)
        

                    

pd_table["Judred_NH2"] = NH2[peptide_numbers].sum(axis=1)
if DoMinMaxScaling:
    pd_table["Judred_NH2"] = (pd_table["Judred_NH2"] / NH2_max) - np.float32(1.0)
        

pd_table["Judred_MW"] = MW[peptide_numbers].sum(axis=1) 
if DoMinMaxScaling:
    pd_table["Judred_MW"] = pd_table["Judred_MW"] - MW_min
    pd_table["Judred_MW"] = pd_table["Judred_MW"] / ((MW_max - MW_min)/2).astype(np.float32)
    pd_table["Judred_MW"] = pd_table["Judred_MW"] - np.float32(1.0)

pd_table["Judred_S"] = S[peptide_numbers].sum(axis=1) 
if DoMinMaxScaling:
    pd_table["Judred_S"] = (pd_table["Judred_S"] / S_max) - np.float32(1.0)

pd_table["Judred_LogP WW"] = (Gwif[peptide_numbers] - Gwoct[peptide_numbers]).sum(axis=1)
if DoMinMaxScaling:
    pd_table["Judred_LogP WW"] = pd_table["Judred_LogP WW"] - LogP_WW_min
    pd_table["Judred_LogP WW"] = pd_table["Judred_LogP WW"] / ((LogP_WW_max - LogP_WW_min)/2.0).astype(np.float32)
    pd_table["Judred_LogP WW"] = pd_table["Judred_LogP WW"] - np.float32(1.0)
        
 
pd_table["Judred_Z"] = charge[peptide_numbers].sum(axis=1) 
if DoMinMaxScaling:
    pd_table["Judred_Z"] = pd_table["Judred_Z"] - Z_min
    pd_table["Judred_Z"] = pd_table["Judred_Z"] / ((Z_max - Z_min)/2.0).astype(np.float32)
    pd_table["Judred_Z"] = pd_table["Judred_Z"] - np.float32(1.0)
        

pd_table["Judred_MaxASA"] = MaxASA[peptide_numbers].sum(axis=1) 
if DoMinMaxScaling:
    pd_table["Judred_MaxASA"] = pd_table["Judred_MaxASA"] - MaxASA_min
    pd_table["Judred_MaxASA"] = pd_table["Judred_MaxASA"] / ((MaxASA_max - MaxASA_min)/2.0).astype(np.float32)
    pd_table["Judred_MaxASA"] = pd_table["Judred_MaxASA"] - np.float32(1.0)

pd_table["Judred_SP2"] = SP2[peptide_numbers].sum(axis=1) 
pd_table["Judred_RotRatio"] = (pd_table["Judred_SP2"]/(SP3[peptide_numbers].sum(axis=1)))
if DoMinMaxScaling:
    pd_table["Judred_RotRatio"] = np.nan_to_num(pd_table["Judred_RotRatio"].values, copy=True)
    pd_table["Judred_RotRatio"] = (pd_table["Judred_RotRatio"] / RotRatio_max) - np.float32(1.0)
    pd_table["Judred_SP2"] = (pd_table["Judred_SP2"] / SP2_max) - np.float32(1.0)
        
pd_table["Judred_Bulkiness"] = bulky[peptide_numbers].sum(axis=1)
if DoMinMaxScaling:
    pd_table["Judred_Bulkiness"] = pd_table["Judred_Bulkiness"] - bulky_min
    pd_table["Judred_Bulkiness"] = pd_table["Judred_Bulkiness"] / ((bulky_max - bulky_min)/2.0).astype(np.float32)
    pd_table["Judred_Bulkiness"] = pd_table["Judred_Bulkiness"] - np.float32(1.0)


pd_table["Judred_OH"] = OH[peptide_numbers].sum(axis=1) 
if DoMinMaxScaling:
    pd_table["Judred_OH"] = (pd_table["Judred_OH"] / OH_max) - np.float32(1.0)

pd_table["Judred_pI"] = pI[peptide_numbers].mean(axis=1)         
if DoMinMaxScaling:
    pd_table["Judred_pI"] = pd_table["Judred_pI"] - pI_min
    pd_table["Judred_pI"] = pd_table["Judred_pI"] / ((pI_max - pI_min)/2.0).astype(np.float32)
    pd_table["Judred_pI"] = pd_table["Judred_pI"] - np.float32(1.0)


#pd_table["Judred_Helical_penalty"] = Helical_penalty[peptide_numbers].sum(axis=1) 
#pd_table["Judred_Helical_penalty"] = (pd_table["Judred_Helical_penalty"] / Helical_penalty_max) - np.float32(1.0)



pd_table.index = peptides
print(pd_table)
pd_table.to_csv("Judred.csv")