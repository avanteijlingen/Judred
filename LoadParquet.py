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
Jparameters["MW"] = Jparameters["MW"].astype(np.float16)
Jparameters["LogP WW"] = Jparameters["LogP WW"].astype(np.float16)
Jparameters["RotRatio"] = Jparameters["RotRatio"].astype(np.float16)
Jparameters["Bulkiness"] = Jparameters["Bulkiness"].astype(np.float16)
print(sys.getsizeof(Jparameters)/1024/1024, "MB")