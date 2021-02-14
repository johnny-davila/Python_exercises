import numpy

from ete3 import NCBITaxa
ncbi = NCBITaxa()

def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 

SP1 = input("Species 1: ")
SP2 = input("Species 2: ")

print(SP1)
print(SP2)

X = ncbi.get_name_translator([SP1])
Y = ncbi.get_name_translator([SP2])

print(X)
print(Y)

X = input("Number 1: ")
Y = input("Number 2: ")

M = ncbi.get_lineage(X)
N = ncbi.get_lineage(Y)
I = intersection(M, N)

A = I[-1]
B = ncbi.get_taxid_translator([A])

print(B)

