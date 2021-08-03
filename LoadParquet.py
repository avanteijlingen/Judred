# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 03:29:47 2021

@author: avtei
"""

import pandas, sys
import numpy as np

Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

L = 6

Jparameters = pandas.read_parquet(Num2Word[L]+"peptides.parquet")

print(Jparameters)

print(sys.getsizeof(Jparameters)/1024/1024, "MB")
Jparameters = Jparameters.astype(np.float16, copy=False)
#Jparameters = Jparameters.astype(np.float32, copy=False)
print(sys.getsizeof(Jparameters)/1024/1024, "MB")

Jdata = Jparameters.values
print(sys.getsizeof(Jdata)/1024/1024, "MB")