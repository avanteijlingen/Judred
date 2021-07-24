# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 03:29:47 2021

@author: avtei
"""

import pandas

Num2Word = {1:"AminoAcids",
            2:"Di",
            3:"Tri",
            4:"Tetra",
            5:"Penta",
            6:"Hexa",
            7:"Hepta",
            8:"Octa"}

L = 3

Jparameters = pandas.read_parquet(Num2Word[L]+"peptides.parquet")

print(Jparameters)