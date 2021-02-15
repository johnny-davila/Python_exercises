from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt 
from Bio import SeqIO
import pandas as pd

print("Your protein multifasta file needs to be in your working directory first.")

IN = input("Write the name of your file here: ")

#Parsing the multifasta protein file and creating a list of ids and sequences
a = SeqIO.parse(IN, "fasta")
ids = []
sequences = []
for i in a:
	ids.append(i.id)
	sequences.append(i.seq)	

#Generating a subset of ids based on the first 1000 sequences of the multifasta
b = ids[:1000]

#Calculating the percentage of aminoacids as specified in the exercise: Cysteine, Acid, Basic, Hydrophobic
a = SeqIO.parse(IN, "fasta")
Cys = []
Aci = []
Bas = []
Hyd = []
for x in a:
	c = ProteinAnalysis(str(x.seq))
	d = float("%0.2f" % c.get_amino_acids_percent()['C'])
	Cys.append(d)
	e = float("%0.2f" % c.get_amino_acids_percent()['D'])
	f = float("%0.2f" % c.get_amino_acids_percent()['E'])
	g = e + f
	Aci.append(g)
	h = float("%0.2f" % c.get_amino_acids_percent()['R'])
	i = float("%0.2f" % c.get_amino_acids_percent()['H'])
	j = float("%0.2f" % c.get_amino_acids_percent()['K'])
	k = g + h + i
	Bas.append(k)
	j = float("%0.2f" % c.get_amino_acids_percent()['K'])
	Hyd.append(j)

#Creating a dataframe based on the values calculated above
DF = pd.DataFrame(list(zip(ids,Cys,Aci,Bas,Hyd)), columns =['Ids','Cysteines','Acid','Basic','Hydrophobic'])

#Generating a dataframe based on the list of ids(b) 
Set = DF[DF["Ids"].isin(b)]

#Creating the boxplot for each type of aminoacid for the entire proteome and the subset of IDs
fig = plt.figure(1, figsize=(9,6))
ax = fig.add_subplot(111)
A = [DF.Cysteines, DF.Acid, DF.Basic, DF.Hydrophobic]
B = [Set.Cysteines, Set.Acid, Set.Basic, Set.Hydrophobic]
bp1 = ax.boxplot(A, positions = [1, 4, 7, 10], widths = 0.6, patch_artist=True, boxprops=dict(facecolor="C0"))
bp2 = ax.boxplot(B, positions = [2, 5, 8, 11], widths = 0.6, patch_artist=True, boxprops=dict(facecolor="C2"))
ax.set_xticklabels(['Cysteines', 'Acid', 'Basic', 'Hydrophobic'])
ax.set_xticks([1.5, 4.5, 7.5, 10.5])
ax.legend([bp1["boxes"][0], bp2["boxes"][0]], ['Total Proteome', 'Subset'], loc='upper right')
fig.savefig('boxcompare.png')

print('The boxplot figure was saved as boxcompare.png in your working directory')
