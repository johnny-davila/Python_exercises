from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D

MW = []
IP = []
II = []

MW2 = []
IP2 = []
II2 = []

markers = []

a = SeqIO.parse("SK_CDS.fasta", "fasta")

for i in a:
	q = ProteinAnalysis(str(i.seq))
	if i.id.find('DNA')>-1:
		r = q.molecular_weight()
		MW.append(r)
		s = q.isoelectric_point()
		IP.append(s)
		t - q.instability_index()
		II.append(t)
	else:
		r = q.molecular_weight()
		MW2.append(r)
		s = q.isoelectric_point()
		IP2.append(s)
		t = q.instability_index()
		II2.append(t) 
	

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(MW, IP, II, c='r', marker=">")
ax.scatter(MW2, IP2, II2, c='b', marker="o")

ax.set_xlabel('Molecular weight')
ax.set_ylabel('Isoelectric point')
ax.set_zlabel('Instability index')

fig.savefig('scatterplot.png')
