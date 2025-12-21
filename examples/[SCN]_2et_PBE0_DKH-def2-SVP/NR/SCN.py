import psi4

psi4.core.set_num_threads(32)
psi4.core.set_output_file('SCN.out', False)
psi4.set_memory('2 GB')

with open('SCN.xyz', 'r') as f:
  xyz_content = f.read().strip()

mol = psi4.geometry(xyz_content)
mol.set_multiplicity(2)
mol.set_molecular_charge(0)

psi4.set_options({
  'basis': 'def2-svp',
  'scf_type': 'df',
  "scf__reference": "uhf",
  'e_convergence': 1e-8,
  'd_convergence': 1e-8,
})

energy, wfn = psi4.energy('pbe0', return_wfn=True)

psi4.molden(wfn, 'SCN.molden')
