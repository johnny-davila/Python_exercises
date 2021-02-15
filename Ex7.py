import numpy
from ete3 import NCBITaxa

#Defining the database of NCBI
ncbi = NCBITaxa()

#Defining a function that will intersect lineages
def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 

#Putting the name of the organisms selected to compare
SP1 = input("Species 1: ")
SP2 = input("Species 2: ")
X = ncbi.get_name_translator([SP1])
Y = ncbi.get_name_translator([SP2])
print(X)
print(Y)

#Putting the code of the organisms selected before
X = input("X: ")
Y = input("Y: ")

#Getting lineages(=ancestors)
M = ncbi.get_lineage(X)
N = ncbi.get_lineage(Y)
I = intersection(M, N)

#Getting the most recent ancestor
A = I[-1]
B = ncbi.get_taxid_translator([A])

print('The most recent common ancestor is:')
print(B)

