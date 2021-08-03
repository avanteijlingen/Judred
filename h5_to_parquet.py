#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 19:14:32 2021

@author: rkb19187
"""
import pandas, h5py, os, sys, time
import peptideutils as pu

L = int(sys.argv[1])
outpath = pu.Num2Word[L].lower()+"peptides.parquet"

if os.path.exists(outpath) and 1 == 0:
    print(outpath, "already exists")
    sys.exit()

h5_file = h5py.File(pu.Num2Word[L].lower()+"peptides.hdf5", 'r')
st = time.time()
names = h5_file["peptides"][()]
print("Loaded peptide names in", round(time.time()-st, 1), "s")
MB_size = sys.getsizeof(names) / 1024 / 1024
print("Names:", MB_size, "MB")
st=time.time()
data = h5_file["data"][()]  #[()]converts them to numpy.ndarray, is stored as float32
print("Loaded data in", round(time.time()-st, 1), "s")
MB_size = sys.getsizeof(data) / 1024 / 1024
print("data:", MB_size, "MB")
st=time.time()
#data = data.astype(np.float32)
Features = h5_file["features"][()]
print("Loaded feature names in", round(time.time()-st, 1), "s")
MB_size = sys.getsizeof(Features) / 1024 / 1024
print("Features:", MB_size, "MB")
h5_file.close()

st=time.time()
Jparameters = pandas.DataFrame(data, index=names, columns=Features)
print("Made dataframe in", round(time.time()-st, 1), "s")
MB_size = sys.getsizeof(Jparameters) / 1024 / 1024
print("Jparameters:", MB_size, "MB")
st=time.time()
Jparameters.to_parquet(outpath)
print("Saved dataframe in", round(time.time()-st, 1), "s")

a="""
h5_file = h5py.File(pu.Num2Word[L].lower()+"peptides.hdf5", 'r')
names = h5_file["peptides"][()]
data = h5_file["data"][()]  #[()]converts them to numpy.ndarray, is stored as float32
#data = data.astype(np.float32)
Features = h5_file["features"][()]
h5_file.close()
Jparameters = pandas.DataFrame(data, index=names, columns=Features)
#parameters = parameters.reindex(targets.index)
Jparameters = Jparameters.sort_index()

#print(Jparameters)


Jparameters.to_parquet(outpath)

print("Saved to:", outpath)
#"""