from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import pandas as pd

print('The multifasta protein file should be in your working directory')

#Parsing the multifasta file
X = input("Write the name of your file here: ")
a = SeqIO.parse(X, "fasta")

#Calculating the percentage of alanines in every sequence
ids = []
seqs = []
Aper = []
Len = []
for i in a:
	b = ProteinAnalysis(str(i.seq))
	c = int(len(str(i.seq)))
	Len.append(c)
	d = "%0.2f" % b.get_amino_acids_percent()['A']
	Aper.append(float(d))
	ids.append(i.id)
	seqs.append(i.seq)

#Calculating the number of sequences corresponding to the 20% of the total sequences
X = int(len(Len))
Y = 0.2*X
Z = int(Y)

#Creating a dataframe and extracting the 20th- and 80th-percentile of the sequences based on alanine content
f = pd.DataFrame(list(zip(ids,seqs,Aper,Len)), columns =['ID','Sequences','Alanine percentage','Length'])
g = f.nlargest(Z, 'Alanine percentage')
h = f.nlargest(X-Z, 'Alanine percentage')

#Generating the two dataframes corresponding to the percentiles
print("The 20th-percentile list is shown below:")
print(g)
print("The 80th-percentile list is shown below:")
print(h)

#Calculating the aminoacids in both percentile groups
T20 = g['Length'].sum()
T80 = h['Length'].sum()
print('The 20th-percentile with higher alanine percentage has', T20, 'aminoacids')
print('The 80th-percentile with higher alanine percentage has', T80, 'aminoacids')
if T20 < T80:
	print('The total of aminoacids is higher in the 80th percentile')
else:
	print('The total of aminoacids is higher in the 20th percentile')
