from RamachanDraw import fetch, phi_psi, plot

#Write the ID of the protein in PDB
X = input("PDB ID: ")

#Generate the Ramachandran for the protein
PDB_id = X
plot(fetch(PDB_id), out='Ram.png')
phi_psi_dict, ignored_res = phi_psi(fetch(PDB_id), return_ignored=True)

print('The Ramachandran plot has been saved as Ram.png')
