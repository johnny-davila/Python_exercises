from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D

print('Your file should be in your working directory and have no spaces in the ID')
X = input("Write the name of your file here: ")

#Parsing the multifasta file
a = SeqIO.parse(X, "fasta")

#Analyzing the aminoacid sequences excluding those containing U or X (as these stop the analysis)
#Separating the analysis for sequences containing the word 'DNA' in their ids
MW = []
IP = []
II = []
MW2 = []
IP2 = []
II2 = []
for i in a:
	if i.seq.find('X')>-1:
		pass
	elif i.seq.find('U')>-1:
		pass
	else:
		if i.id.find('DNA')>-1:
			q = ProteinAnalysis(str(i.seq))
			r = q.molecular_weight()
			MW.append(r)
			s = q.isoelectric_point()
			IP.append(s)
			t = q.instability_index()
			II.append(t)
		else:
			m = ProteinAnalysis(str(i.seq))
			u = m.molecular_weight()
			MW2.append(u)
			v = m.isoelectric_point()
			IP2.append(v)
			w = m.instability_index()
			II2.append(w)

#Creating the 3D scatterplot based on those two sets of information
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.text2D(0.05, 0.95,"DNA (red) vs. non-DNA (blue) proteins", transform=ax.transAxes)
ax.scatter(MW, IP, II, c='r', marker=">")
ax.scatter(MW2, IP2, II2, c='b', marker="o")
ax.set_xlabel('Molecular weight')
ax.set_ylabel('Isoelectric point')
ax.set_zlabel('Instability index')
fig.savefig('3Dscatterplot.png')
print('The figure was saved as 3Dscatterplot.png')
