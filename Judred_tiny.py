# -*- coding: utf-8 -*-
"""

@author: avtei
"""


import pandas, sys, os, time, math
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from datetime import datetime

try:
    import cupy as cp
    mempool = cp.get_default_memory_pool()
    mempool.set_limit(size=3*1024**3) # 3 GB
    use_gpu = True
except:
    print("Warning: not able to use GPU, is cupy installed?")
    use_gpu = False


Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

L = int(sys.argv[1])
#chunksize = min([math.floor((20**L)/2), 2560000])
chunksize = int(sys.argv[2])
if len(sys.argv) > 3:
    assert sys.argv[3].lower() == "peptoids"
    peptoids = True
else:
    peptoids = False
    

if peptoids:
    peptoid_data = pandas.read_csv("peptoid_data/Peptoid_Data.csv", index_col=0)
    fname = Num2Word[L].lower()+"peptoids_normalized.parquet"
    letters_1 = peptoid_data.index
    SP2 = peptoid_data["SP2"].values
    SP3 = peptoid_data["SP3"].values
    NH2 = peptoid_data["NH2"].values
    MW = peptoid_data["Mass (AMU)"].values
    charge = peptoid_data["Charge"].values
    S = peptoid_data["S"].values
    Gwif = peptoid_data["Gwif"].values
    Gwoct = peptoid_data["Gwoct"].values
    MaxASA = peptoid_data["ASA"].values
    OH = peptoid_data["OH"].values
    features = ["SP2", "NH2", "MW", "S", "Solubility Scale", "Z", "MaxASA", "RotRatio", "OH"]
    
    #Check the chunksize is divisible
    assert (peptoid_data.shape[0]**L) % chunksize == 0, "peptoid_data.shape[0]**L is not devisible by the chunksize"
else:
    fname = Num2Word[L].lower()+"peptides_normalized.parquet"
    letters_1 = np.array(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
    
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
numbers = np.arange(0, len(letters_1), dtype=np.uint8)
steps = np.array([numbers.shape[0]**i for i in range(L-1, -1, -1)], dtype=np.uint64)

if use_gpu:
    steps_gpu = cp.array(steps)
    SP2_gpu = cp.array(SP2)
    SP3_gpu = cp.array(SP3)
    NH2_gpu = cp.array(NH2)
    MW_gpu = cp.array(MW)
    S_gpu = cp.array(S)
    charge_gpu = cp.array(charge)
    Gwif_gpu = cp.array(Gwif)
    Gwoct_gpu = cp.array(Gwoct)
    MaxASA_gpu = cp.array(MaxASA)
    OH_gpu = cp.array(OH)
    if not peptoids:
        pI_gpu = cp.array(pI)
        bulky_gpu = cp.array(bulky)

if L <= 4:
    indices = np.indices((len(numbers),) * L, dtype=np.uint8, sparse=False)
    indices = np.rollaxis(indices, 0, L + 1)
    indices = indices.reshape(-1, L)


SP2_max = ((max(SP2)*L)/2.0).astype(np.float32) 
if peptoids:
    polytryptophan_index = [18]*L
    RotRatio_max = SP2[polytryptophan_index].sum() / SP3[polytryptophan_index].sum() 
    RotRatio_max = np.float32(RotRatio_max/2.0)
else:
    RotRatio_max = 10.0*L
    
NH2_max = ((max(NH2)*L)/2.0).astype(np.float32)
MW_min = (min(MW)*L).astype(np.float32)
MW_max = (max(MW)*L).astype(np.float32)
S_max = ((max(S)*L)/2.0).astype(np.float32)
Z_min = (min(charge)*L).astype(np.float32)
Z_max = (max(charge)*L).astype(np.float32)
if not peptoids:
    polyasparticacid_index = [2]*L
    polyisoleucine_index = [7]*L
    bulky_min = (min(bulky)*L).astype(np.float32)
    bulky_max = (max(bulky)*L).astype(np.float32)
    pI_min = (min(pI)).astype(np.float32)
    pI_max = (max(pI)).astype(np.float32)
else:
    polyasparticacid_index = [np.where(letters_1 == "Nd")[0][0]]*L
    polyisoleucine_index = [np.where(letters_1 == "Ni")[0][0]]*L


LogP_WW_min = (Gwif[polyasparticacid_index] - Gwoct[polyasparticacid_index]).sum()
LogP_WW_max = (Gwif[polyisoleucine_index] - Gwoct[polyisoleucine_index]).sum()
MaxASA_min = (min(MaxASA)*L).astype(np.float32)
MaxASA_max = (max(MaxASA)*L).astype(np.float32)
OH_max = ((max(OH)*L)/2.0).astype(np.float32)


if use_gpu:
    SP2_gpu_max = ((max(SP2_gpu)*L)/2.0).astype(np.float32) 
    RotRatio_gpu_max = SP2_gpu[polytryptophan_index].sum() / SP3_gpu[polytryptophan_index].sum() 
    RotRatio_gpu_max = RotRatio_gpu_max/2.0
    NH2_gpu_max = ((max(NH2_gpu)*L)/2.0)
    MW_gpu_min = (min(MW_gpu)*L).astype(cp.float32)
    MW_gpu_max = (max(MW_gpu)*L).astype(cp.float32)
    S_gpu_max = ((max(S_gpu)*L)/2.0).astype(cp.float32)
    Z_gpu_min = (min(charge_gpu)*L).astype(cp.float32)
    Z_gpu_max = (max(charge_gpu)*L).astype(cp.float32)
    LogP_WW_gpu_min = (Gwif_gpu[polyasparticacid_index] - Gwoct_gpu[polyasparticacid_index]).sum()
    LogP_WW_gpu_max = (Gwif_gpu[polyisoleucine_index] - Gwoct_gpu[polyisoleucine_index]).sum()
    MaxASA_gpu_min = (min(MaxASA_gpu)*L).astype(cp.float32)
    MaxASA_gpu_max = (max(MaxASA_gpu)*L).astype(cp.float32)
    OH_gpu_max = ((max(OH_gpu)*L)/2.0).astype(cp.float32)
    if not peptoids:
        bulky_gpu_min = (min(bulky_gpu)*L).astype(cp.float32)
        bulky_gpu_max = (max(bulky_gpu)*L).astype(cp.float32)
        pI_gpu_min = (min(pI)).astype(cp.float32)
        pI_gpu_max = (max(pI)).astype(cp.float32)



y = np.zeros((chunksize, len(features)), dtype=np.float32)
pd_table = pandas.DataFrame(y, columns=features, index=np.arange(0,chunksize))
print(pd_table)
pd_table["SP2"] = pd_table["SP2"].astype(np.float32)
pd_table["NH2"] = pd_table["NH2"].astype(np.float32)
pd_table["MW"] = pd_table["MW"].astype(np.float32)
pd_table["S"] = pd_table["S"].astype(np.float32)
pd_table["Z"] = pd_table["Z"].astype(np.float32)
pd_table["MaxASA"] = pd_table["MaxASA"].astype(np.float32)
pd_table["RotRatio"] = pd_table["RotRatio"].astype(np.float32)
pd_table["OH"] = pd_table["OH"].astype(np.float32)
if not peptoids:
    pd_table["pI"] = pd_table["pI"].astype(np.float32)
    pd_table["Bulkiness"] = pd_table["Bulkiness"].astype(np.float32)
else:
    SolubilityColumnName = "Solubility Scale"
pd_table[SolubilityColumnName] = pd_table[SolubilityColumnName].astype(np.float32)


table = pa.Table.from_pandas(pd_table, preserve_index=False)

# =============================================================================
# print(table.schema)
# sys.exit()
# =============================================================================

times = []
#use_gpu = False
print("Use GPU:", use_gpu)
with pq.ParquetWriter(fname, table.schema) as writer:
    iterations = (numbers.shape[0]**L)//chunksize
    iterations += int(bool((numbers.shape[0]**L)%chunksize))
    for iteration in range(iterations):
        index_lower = iteration*chunksize
        index_upper = (iteration+1)*chunksize
        if index_upper > numbers.shape[0]**L:
            chunksize = numbers.shape[0]**L - index_lower
            index_upper = numbers.shape[0]**L
            y = np.zeros((chunksize, len(features)), dtype=np.float32)
        print(index_lower, "-", index_upper, "|", numbers.shape[0]**L, chunksize)
        #print("="*45)
        #continue
        st = time.time()
        # Where the min value is 0 we can do a faster normalization

        pd_table = pandas.DataFrame(y, columns=features, index=np.arange(0,chunksize))
        
        # For the larger datasets the indices become to large so we have to chunk that too
        if L > 4:
            if use_gpu:
                index = cp.arange(index_lower, index_upper, dtype=np.uint64)
                peptide_numbers_gpu = cp.zeros((chunksize, L), dtype=np.uint8)
                for i in range(L):
                    peptide_numbers_gpu[:,i] = index // steps[i]
                    v = (index//steps_gpu[i])
                    v = v * steps_gpu[i]
                    index = index - v
            else:
                index = np.arange(index_lower, index_upper, dtype=np.uint64)
                peptide_numbers = np.zeros((chunksize, L), dtype=np.uint8)
                for i in range(L):
                    peptide_numbers[:,i] = index // steps[i]
                    v = (index//steps[i])
                    v = v * steps[i]
                    index = index - v
        else:
            if use_gpu:
                peptide_numbers_gpu = cp.array(indices[index_lower:index_upper])
            else:
                peptide_numbers = indices[index_lower:index_upper]
                
        if use_gpu:
            chunk = NH2_gpu[peptide_numbers_gpu]
            chunk = chunk.sum(axis=1)
            chunk = (chunk / NH2_gpu_max ) - 1
            pd_table["NH2"] = chunk.get()
        else:
            pd_table["NH2"] = NH2[peptide_numbers].sum(axis=1)
            pd_table["NH2"] = (pd_table["NH2"] / NH2_max) - np.float32(1.0)
        
        if use_gpu:
            chunk = MW_gpu[peptide_numbers_gpu].sum(axis=1) 
            chunk = chunk - MW_gpu_min
            chunk = chunk / ((MW_gpu_max - MW_gpu_min)/2)
            chunk = chunk - 1
            pd_table["MW"] = chunk.get()
        else:
            pd_table["MW"] = MW[peptide_numbers].sum(axis=1) 
            pd_table["MW"] = pd_table["MW"] - MW_min
            pd_table["MW"] = pd_table["MW"] / ((MW_max - MW_min)/2).astype(np.float32)
            pd_table["MW"] = pd_table["MW"] - np.float32(1.0)
        
        if use_gpu:
            chunk = S_gpu[peptide_numbers_gpu]
            chunk = chunk.sum(axis=1)
            chunk = (chunk / S_gpu_max ) - 1
            pd_table["S"] = chunk.get()
        else:
            pd_table["S"] = S[peptide_numbers].sum(axis=1) 
            pd_table["S"] = (pd_table["S"] / S_max) - np.float32(1.0)
        
        if use_gpu:
            logp_gpu_gwif = Gwif_gpu[peptide_numbers_gpu].sum(axis=1)
            chunk = Gwoct_gpu[peptide_numbers_gpu].sum(axis=1)
            chunk = logp_gpu_gwif - chunk
            chunk = chunk - LogP_WW_min
            chunk = chunk / ((LogP_WW_max - LogP_WW_min)/2.0)
            chunk = chunk - 1
            pd_table[SolubilityColumnName] = chunk.get()
        else:
            pd_table[SolubilityColumnName] = (Gwif[peptide_numbers] - Gwoct[peptide_numbers]).sum(axis=1)
            pd_table[SolubilityColumnName] = pd_table[SolubilityColumnName] - LogP_WW_min
            pd_table[SolubilityColumnName] = pd_table[SolubilityColumnName] / ((LogP_WW_max - LogP_WW_min)/2.0).astype(np.float32)
            pd_table[SolubilityColumnName] = pd_table[SolubilityColumnName] - np.float32(1.0)
        
        if use_gpu:
            chunk = charge_gpu[peptide_numbers_gpu]
            chunk = chunk.sum(axis=1)
            chunk = chunk - Z_gpu_min
            chunk = chunk / ((Z_max - Z_min)/2.0)
            chunk = chunk - 1
            pd_table["Z"] = chunk.get()
        else:
            pd_table["Z"] = charge[peptide_numbers].sum(axis=1) 
            pd_table["Z"] = pd_table["Z"] - Z_min
            pd_table["Z"] = pd_table["Z"] / ((Z_max - Z_min)/2.0).astype(np.float32)
            pd_table["Z"] = pd_table["Z"] - np.float32(1.0)
        
        if use_gpu:
            chunk = MaxASA_gpu[peptide_numbers_gpu]
            chunk = chunk.sum(axis=1)
            chunk = chunk - MaxASA_gpu_min
            chunk = chunk / ((MaxASA_gpu_max - MaxASA_gpu_min) / 2.0)
            chunk = chunk - 1
            pd_table["MaxASA"] = chunk.get()
        else:
            pd_table["MaxASA"] = MaxASA[peptide_numbers].sum(axis=1) 
            pd_table["MaxASA"] = pd_table["MaxASA"] - MaxASA_min
            pd_table["MaxASA"] = pd_table["MaxASA"] / ((MaxASA_max - MaxASA_min)/2.0).astype(np.float32)
            pd_table["MaxASA"] = pd_table["MaxASA"] - np.float32(1.0)
        
        if use_gpu:
            SP2_chunk = SP2_gpu[peptide_numbers_gpu].sum(axis=1)
            chunk = SP3_gpu[peptide_numbers_gpu].sum(axis=1)
            chunk = SP2_chunk / chunk
            #print((cp.isnan(RotRatio_chunk) == True).any())
            chunk = cp.nan_to_num(chunk)
            chunk = (chunk / RotRatio_gpu_max) - 1
            SP2_chunk = (SP2_chunk / SP2_gpu_max) - 1
            pd_table["SP2"] = SP2_chunk.get()
            pd_table["RotRatio"] = chunk.get()
        else:
            pd_table["SP2"] = SP2[peptide_numbers].sum(axis=1) 
            pd_table["RotRatio"] = (pd_table["SP2"]/(SP3[peptide_numbers].sum(axis=1)))
            pd_table["RotRatio"] = np.nan_to_num(pd_table["RotRatio"].values, copy=True)
            pd_table["RotRatio"] = (pd_table["RotRatio"] / RotRatio_max) - np.float32(1.0)
            pd_table["SP2"] = (pd_table["SP2"] / SP2_max) - np.float32(1.0)
        
        if use_gpu:
            chunk = OH_gpu[peptide_numbers_gpu]
            chunk = chunk.sum(axis=1)
            chunk = (chunk / OH_gpu_max ) - 1
            pd_table["OH"] = chunk.get()
        else:
            pd_table["OH"] = OH[peptide_numbers].sum(axis=1) 
            pd_table["OH"] = (pd_table["OH"] / OH_max) - np.float32(1.0)
            
        if not peptoids:
            if use_gpu:
                chunk = bulky_gpu[peptide_numbers_gpu]
                chunk = chunk.sum(axis=1)
                chunk = chunk - bulky_gpu_min
                chunk = chunk / ((bulky_gpu_max - bulky_gpu_min)/2.0)
                chunk = chunk - 1
                pd_table["Bulkiness"] = chunk.get()
            else:
                pd_table["Bulkiness"] = bulky[peptide_numbers].sum(axis=1)
                pd_table["Bulkiness"] = pd_table["Bulkiness"] - bulky_min
                pd_table["Bulkiness"] = pd_table["Bulkiness"] / ((bulky_max - bulky_min)/2.0).astype(np.float32)
                pd_table["Bulkiness"] = pd_table["Bulkiness"] - np.float32(1.0)
            if use_gpu:
                chunk = pI_gpu[peptide_numbers_gpu]
                chunk = chunk.mean(axis=1)
                chunk = (chunk / pI_gpu_max ) - 1
                pd_table["pI"] = chunk.get()
            else:
                pd_table["pI"] = pI[peptide_numbers].mean(axis=1) 
                pd_table["pI"] = (pd_table["pI"] / pI_max) - np.float32(1.0)

        
        table = pa.Table.from_pandas(pd_table.astype(np.float32), preserve_index=False)
        writer.write_table(table)
        
        times.append(time.time()-st)
        eta = sum(times) / (iteration+1)
        eta = (eta*iterations) - (eta*(iteration+1))
        eta_str = datetime.fromtimestamp(eta).strftime("%H:%M:%S")
        print("{:.2e}".format(index_upper), "/", "{:.2e}".format(numbers.shape[0]**L), round(times[-1], 3), "s eta:", eta_str)
        
        #break

if use_gpu:
    del peptide_numbers_gpu
    del SP2_chunk
    del logp_gpu_gwif
    del chunk
    cp._default_memory_pool.free_all_blocks()



print(Num2Word[L].lower()+"peptides done in", round(sum(times), 3), "s use_gpu:", use_gpu)