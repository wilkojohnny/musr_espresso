# qe_parameterfinder.py - create and run quantum espresso files, and from them produce plottable files to determine
# the optimum parameters to use
# John Wilkinson 6/5/19
#
# a_sweep(..): sweep up the lattice parameter, and observe how the total energy changes
# sweep(..): sweep up any other parameter in the pw.x input file, outputting parameter and energy
# get_energy(..): calculate the total energy using pw.x

from ase.build import bulk  # to set up the cells
from ase import io, calculators  # to get the structures from CIF files, and for error handling
from ase.calculators.espresso import Espresso  # to do espresso calculations
import gle_utils  # for plotting with GLE
import numpy as np  # so we can do arange (i.e same as range() but with floats)
import copy  # to create copies of the structure


# define pseudopotential file names
pseudopotentials = { 'C': 'C.pbe-n-rrkjus_psl.1.0.0.UPF',
                     'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
                     'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
                     'F': 'F.pbe-n-rrkjus_psl.0.1.UPF'}


def main():

    input_data = {
        'system': {
             'ecutwfc': 85,
             'ecutrho': 750,
             #'lda_plus_u': True,
             #'Hubbard_U(1)': 3,
             # 'occupations': 'smearing',
             # 'smearing': 'gaussian',
             # 'degauss': 0.01
                },
        'electrons': {
            'mixing_beta': 0.7,
            'mixing_mode': 'plain',
        }
    }

    cif_location = 'PERFECTA.cif'

    # get the atoms
    atoms = io.read(cif_location)
    close_atoms = copy.deepcopy(atoms)
    close_atoms.set_cell(atoms.get_cell()*0.95, scale_atoms=True)

    #lattice_parameter_sweep(0.8, 1.2, 0.02, atoms, input_data, 3)

    all_parameters = []
    dE = []
    for nk in range(2, 6):
        parameters = np.arange(40, 75, 5)
        parameters, energies = sweep(atoms, input_data, 'system', 'ecutwfc', nk, parameters)
        close_atoms_parameters, close_energy = sweep(close_atoms, input_data, 'system', 'ecutwfc', nk, parameters)
        all_parameters.append(close_atoms_parameters)

        this_dE = []
        # subtract the unperturbed lattice energies from the perturbed (small) - but only subtract if not none...
        # for each close_atom_parameter
        for i_clp in range(0, len(close_atoms_parameters)):
            clp = close_atoms_parameters[i_clp]
            # find the location in the array of parameters where this cl_par is
            i_p = parameters.index(clp)
            # use that index to find the energy of the normal-sized (unperturbed) lattice at this parameter
            energy_change = energies[i_p] - close_energy[i_clp]
            # append to dE
            this_dE.append(energy_change)

        dE.append(this_dE)

    gle_utils.plot_xy(all_parameters, dE, 'ecutwfc and scf energy, PERFECTA', legend=["nk=" + str(i) for i in range(2, 6)])

    return 0


# sweep up the lattice parameters, and plot a graph
# sweeps up from scale_min to scale_max in steps of scale_step
# postamble of "! plot the experimental lattice parameter \n "
#             + "set color red \n amove xg(4.62) yg(ygmin) \n aline xg(4.62) yg(ygmax)"
# is a nice way to plot the experimental value of the lattice parameter
def lattice_parameter_sweep(scale_min, scale_max, scale_step, atoms, input_data, nk, postamble='', plot=True):
    a = []
    e = []

    original_cell = atoms.get_cell()

    current_atoms = copy.deepcopy(atoms)

    for scale in np.arange(scale_min, scale_max, scale_step):
        # check scale is not 0
        assert scale != 0

        # alter
        current_atoms.set_cell(scale*original_cell, scale_atoms=True)

        energy = get_energy(atoms=current_atoms, nk=nk, input_data=input_data)/13.6056980659  # gets converted into Ry

        a.append(scale)
        e.append(energy)

        print(str(scale) + ' ' + str(energy))


    # plot
    if plot:
        gle_utils.plot_xy(a, e, title="DFT total energy vs Lattice scale", legend="DFT Energy /Ry",
                          postamble="! plot the experimental lattice parameter \n "
                                    + "set color red \n amove xg(1) yg(ygmin) \n aline xg(1) yg(ygmax)" + postamble)

    return a, e


def sweep(atoms, input_data, sweep_namespace, sweep_parameter, nk, parameters):
    successful_parameters = []
    e = []

    # for each parameter value...
    for parameter_value in parameters:
        # ... set up the input data
        try:
            # try to just write the new value into the dictionary - if it doesn't exist, add it
            input_data[sweep_namespace][sweep_parameter] = parameter_value
        except KeyError:
            # add the key and the value if not already in it
            input_data[sweep_namespace] = {sweep_parameter: parameter_value}

        # calculate the energy and add to the corresponding array
        energy = get_energy(atoms, nk, input_data)

        # if energy is not none (i.e the QE run converged)
        if energy is not None:
            e.append(energy)
            # save the current parameter value to the array
            successful_parameters.append(parameter_value)
            # print out the parameters and the calculated energy
            print(str(nk) + ' ' + str(parameter_value) + ' ' + str(energy))
        else:
            # energy is none: warn the user
            print('failed: ' + str(nk) + ' ' + str(parameter_value) + ' ' + str(energy))

    return successful_parameters, e


def get_energy(atoms: bulk, nk=3, input_data=None):
    # fix mutable arguments
    if input_data is None:
        input_data = {}

    # set up the calculator
    calc = Espresso(pseudopotentials=pseudopotentials,
                    tstress=False, tprnfor=False, kpts=(nk, nk, nk), input_data=input_data)
    # attach the calculator to the atoms
    atoms.set_calculator(calc)

    # calculate energy with DFT, and return energy
    try:
        energy = atoms.get_total_energy()
    except calculators.calculator.CalculationFailed:
        energy = None

    return energy

if __name__ == '__main__':
    main()
