from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

a = SeqIO.parse("SK_CDS.fasta", "fasta")

ids = []
Aper = []
Len = []

for i in a:
	b = ProteinAnalysis(str(i.seq))
	c = int(len(str(i.seq)))
	Len.append(c)
	d = "%0.2f" % b.get_amino_acids_percent()['A']
	Aper.append(d)
	ids.append(i.id)

X = int(len(Len))
Y = 0.2*X
Z = int(Y)

f = pd.DataFrame(list(zip(ids,Len)), columns =['ID','Length'])
g = f.nlargest(Z, 'Length')
h = f.nlargest(X-Z, 'Length')

T20 = g['Length'].sum()
T80 = h['Length'].sum()

print(T20)
print(T80)

if T20 < T80:
	print('The total of aminoacids is higher in the 80th percentile')
else:
	print('The total of aminoacids is higher in the 20th percentile')
