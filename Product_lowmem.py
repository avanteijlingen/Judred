# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 14:07:05 2021

@author: rudolf
"""
import numpy as np
import itertools, random, time

L = 4
steps = np.array([20**i for i in range(L-1, -1, -1)], dtype=np.uint64)
print(steps)
chunksize = 10

numbers = np.arange(0, 20)
number_arrs = [numbers]*L


Validation = np.array(list(itertools.product(*number_arrs)))


st = time.time()
index = np.arange(0, 20**L)
solution = np.zeros((20**L, L), dtype=np.uint8)
for i in range(L):
    solution[:,i] = index // steps[i]
    v = (index//steps[i])
    v = v * steps[i]
    index = index - v
tri = time.time()-st
penta = round(tri * 20 * 20, 2)


print("penta:", penta, "s")

for i, j in zip(Validation, solution):
    print(i, j, (i==j).all())
    