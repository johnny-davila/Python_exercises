from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt 
from Bio import SeqIO
import pandas as pd

a = SeqIO.parse("Prot_ZF.fasta", "fasta")

ids = []
sequences = []
for i in a:
	ids.append(i.id)
	sequences.append(i.seq)	
	
b = ids[:1000]

a = SeqIO.parse("Prot_ZF.fasta", "fasta")

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

DF = pd.DataFrame(list(zip(ids,Cys,Aci,Bas,Hyd)), columns =['Ids','Cysteines','Acid','Basic','Hydrophobic'])

Set = DF[DF["Ids"].isin(b)]

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

#Cysteines
#Acid
#Basic
#Hydrophobic
