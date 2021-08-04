#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 16:47:52 2021

@author: rkb19187
"""
import numpy as np
import sys, itertools


L = 2

def pep2index(peptide):
    L = len(peptide)
    size = int(20**L)
    solution = 0
    letters_1 = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    for i in range(1, L+1):
        index = np.where(letters_1 == peptide[i-1])[0][0]
        number = int((size/(20**i)) * index)
        #print(index, number)
        solution += number
    return solution
        
        
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

letters_1 = list("ACDEFGHIKLMNPQRSTVWY")
letters_set = [letters_1]*L
Validation = ["".join(x) for x in list(itertools.product(*letters_set))]
a="""
for i in [2105]:
    print(i)
    print(Validation[i], end=" - ")
    pep = index2pep(i, L)
    print(pep)
    index = pep2index(pep)
    print(index)
#"""
    
for pep in Validation:
    index = pep2index(pep)
    print(pep, index)
