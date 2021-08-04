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

L = 3

Jparameters = pandas.read_parquet(Num2Word[L].lower()+"peptides_NORMALIZED.parquet")
#
print(Jparameters)
if Jparameters.shape[0] == 0:
    sys.exit()
print(sys.getsizeof(Jparameters)/1024/1024, "MB")

for col in Jparameters.columns:
    print(col, Jparameters[col].max(), Jparameters[col].min())

    
a="""
Jparameters = Jparameters.astype(np.float16)
print(sys.getsizeof(Jparameters)/1024/1024, "MB")
"""