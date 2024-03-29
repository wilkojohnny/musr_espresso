"""
NaF_ecutwfc.py -- example for using musr_espresso's PW() class to find the corrrect cutoff parameters for NaF
"""

import musr_espresso

# import the structure from the CIF file
NaF_atoms = musr_espresso.io.read('NaF.cif')

# give the pseudopotential file names
pseudopotentials = { 'Na': 'Na.pbe-spnl-rrkjus_psl.1.0.0.UPF',
                     'F': 'F.pbe-n-rrkjus_psl.0.1.UPF'}

# give the pw.x parameters
input_data = {
    'control': {
        'pseudo_dir': 'pp',
        'outdir': 'out'
    },
    'system': {
        'ecutwfc': 85,
        'ecutrho': 900,
    },
    'electrons': {
        'mixing_beta': 0.7,
    }
}

# create an instance of the PW class
pw = musr_espresso.PW(pwi_params=input_data, atoms=NaF_atoms, pseudopotentials=pseudopotentials, nk=(3,3,3))

# scan up ecutwfc and nk
nks = [(nk, nk, nk) for nk in range(2, 5)]
param_energies = pw.sweep_parameter('system', 'ecutwfc', range(30, 75, 5), nk=nks, plot=True, small_sf=(0.8, 1, 1))

# print the result
print(param_energies)