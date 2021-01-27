"""
NaF_bands.py -- use musr_espresso to calculate the band structure of NaF
"""

import musr_espresso
import numpy as np

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

# define the points on the k-path:
k_points = [["W", np.array((1/4, 3/4, 1/2))],
            ["L", np.array((1/2, 1/2, 1/2))],
            ["\\Gamma", np.array((0, 0, 0))],  # to do LaTEX characters, do \\ instead of just \
            ["X", np.array((0, 1/2, 1/2))],
            ["W", np.array((1/4, 3/4, 1/2))],
            ["K", np.array((3/8, 3/4, 3/8))]]

# run the calculation!
pw.calculate_band_structure(k_points, 10, 40, cleanup=True)
