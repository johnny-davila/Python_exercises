from RamachanDraw import fetch, phi_psi, plot

PDB_id = '1MBN'

plot(fetch(PDB_id))

phi_psi_dict, ignored_res = phi_psi(fetch(PDB_id), return_ignored=True)
