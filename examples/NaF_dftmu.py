"""
NaF_dftmu.py -- example for using musr_espresso's PW() class to set up DFT+mu runs, ready for the cluster
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

pw.setup_dftmu(muon_positions=10, supercell=[2,2,2], run_time=12*60*60, nodes=3, output_directory='NaF',
               prefix='NaF.relax.mu')