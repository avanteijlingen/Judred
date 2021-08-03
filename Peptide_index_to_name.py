#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 16:47:52 2021

@author: rkb19187
"""
import peptideutils as pu
import numpy as np
import sys

L = 4
size = 20**L
step = int(size/20)
above_size = int(size/20)
steps = np.array([i*step for i in range(20)])
Ls = pu.peptideutils_letters1

peptides = pu.GenerateDatasetIndex(L)

test = 999
if test >= size:
    print("No peptide this large in the dataset")
    sys.exit()
    
print(peptides[test], end=" - ")

solution = []

letter_i = np.where(steps <= test)[0][-1]
letter = pu.peptideutils_letters1[letter_i]
solution.append(letter)

while len(solution) < L:
    test = test%step
    step = int(step/20)
    steps = np.array([i*step for i in range(20)])
    if test <= above_size/20 and len(solution) < L-1:
        letter = "A"
        #print(test <= above_size/20)
    else:
        letter_i = np.where(steps <= test)[0][-1]
        letter = pu.peptideutils_letters1[letter_i]
        #print(letter_i, letter)
    above_size = int(above_size/20)
    
    solution.append(letter)


print("".join(solution))