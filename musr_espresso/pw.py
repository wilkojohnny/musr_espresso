"""
pw.py -- unit to do pw.x calculations using a nice, clean, simple python interface

Created by John Wilkinson, 8/10/2020
"""

from ase import io, build, atoms
from ase.visualize import view
from ase.build import bulk
from ase.calculators import espresso, calculator
import copy
import numpy as np
import sys
import time
from .pointgrid import point_grid
from .gle_utils import plot_xy, plot_scatter, plot_bands
import subprocess


import functools
print = functools.partial(print, flush=True)


class PW(object):
    """
    PW class -- contains everything you could ever want to do with Quantum Espresso's pw.x utility
    """

    def __init__(self, pwi_params: dict, atoms: bulk, pseudopotentials: dict, nk: tuple):
        print('☕️ Welcome to musr_espresso -- DFT+mu tools for Quantum Espresso\n')
        mpi_args = ' '.join(sys.argv[1:])
        if mpi_args == '':
            print('😴 Using serial mode. To run in parallel, just add -n <n_processors> to the end of the python '
                  'command')
        else:
            print('😁 Using parallel mode. MPI options given are ' + mpi_args + '\n')
        self.pwi_params = pwi_params
        self.atoms = atoms
        self.pseudopotentials = pseudopotentials
        if isinstance(nk, tuple):
            self.nk = nk
        else:
            self.nk = (nk, nk, nk)

    def sweep_lattice_parameter(self, start_scale, end_scale, scale_step, plot=False):
        """
        sweep the lattice parameter from (a,b,c)*start to (a,b,c)*stop, going in steps of (a,b,c)*step, using pw.x to
        calculate the energy for each
        :param start_scale: starting lattice parameter/(a,b,c)
        :param end_scale: final lattice parameter/(a,b,c)
        :param scale_step: incrementation of the lattice parameter
        :param plot: do a plot of the data
        :return: 2D array, [0] = lattice parameter, [1] = energy
        """
        original_cell = self.atoms.get_cell()

        current_atoms = copy.deepcopy(self.atoms)

        a, e = [], []

        for scale in np.arange(start_scale, end_scale, scale_step):
            # check scale is not 0
            assert scale != 0

            # alter the lattice parameter of the cell
            current_atoms.set_cell(scale * original_cell, scale_atoms=True)

            energy = self.get_energy(atoms=current_atoms) / 13.6056980659  # convert into rydbergs

            a.append(scale)
            e.append(energy)

            print(str(scale) + ' ' + str(energy))

        # plot
        if plot:
            plot_xy(a, e, title="DFT total energy vs Lattice scale", legend="DFT Energy /Ry",
                    postamble="! plot the experimental lattice parameter \n "
                               + "set color red \n amove xg(1) yg(ygmin) \n aline xg(1) yg(ygmax)",
                    xtitle="lattice scale factor", ytitle="DFT Energy (Ry)")
        return [a, e]

    def sweep_parameter(self, pwi_namespace, pwi_parameter, pwi_param_values, nk=None, small_sf=0.95, plot=False):
        """
        Sweep a parameter in the pwi file, calculating E(a)-E(small_sf*a) for each
        :param pwi_namespace: namespace the parameter to change is in
        :param pwi_parameter: parameter to change
        :param pwi_param_values: array of values to try
        :param nk: array of nks to try (if None, just uses the class's saved nk (i.e self.nk))
        :param small_sf: scale factor to apply to the cell for the small cell energy calculation
        :param plot: do plot of results
        :return: list with entries [nk_id][param][dE]
        """

        # if nk is a list, then loop up it
        if isinstance(nk, list):
            result = []
            for i_nk, this_nk in enumerate(nk):
                # get one of the list elements, and re-call this function with *one* of the k-points, not plotting
                result.append(self.sweep_parameter(pwi_namespace, pwi_parameter, pwi_param_values, this_nk, plot=False))
            if plot:
                all_parameters = [res[0] for res in result]
                dE = [res[1] for res in result]
                plot_xy(all_parameters, dE, pwi_parameter + ' and scf energy,' + str(self.atoms.symbols),
                        legend=["nk=" + str(i) for i in nk], ytitle="E(a)-E({:.3f}a) (Ry)".format(small_sf),
                        xtitle=pwi_parameter)
            return result
        elif nk is None:
            # nk's default is self.nk
            nk = self.nk
        else:
            # if nk is just one integer, assume the user wants a k-grid of (nk, nk, nk)
            if isinstance(nk, int):
                nk = (nk, nk, nk)
            parameters, energies = self.do_sweep(sweep_param_namespace=pwi_namespace,
                                                 sweep_parameter=pwi_parameter,
                                                 nk=nk,
                                                 parameters=pwi_param_values)

            close_atoms = copy.deepcopy(self.atoms)
            close_atoms.set_cell(self.atoms.get_cell()*small_sf, scale_atoms=True)

            close_atoms_parameters, close_energy = self.do_sweep(atoms=close_atoms,
                                                                 sweep_param_namespace=pwi_namespace,
                                                                 sweep_parameter=pwi_parameter,
                                                                 nk=nk,
                                                                 parameters=parameters)

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

            print('nk=' + str(nk))
            print(this_dE)

            if plot:
                plot_xy(close_atoms_parameters, this_dE, pwi_parameter + ' and scf energy,' +
                        str(self.atoms.symbols), legend="k-grid=" + str(nk),
                        ytitle="E(a)-E({:.3f}a) (Ry)".format(small_sf), xtitle=pwi_parameter)

            return [close_atoms_parameters, this_dE]

    def get_energy(self, atoms: bulk = None, nk=None, pwi_params=None):
        """
        get_energy() : use pw.x to calculate the energy of the structure
        :param atoms: ASE atoms to calculate the energy of
        :param nk: tuple of k points (nk, nk, nk)
        :param pwi_params: parameters in the pwi file to use (uses default if not specified)
        :return: energy in Ry
        """
        # if anything is None, use the class's definition
        if atoms is None:
            atoms = self.atoms
        if nk is None:
            nk = self.nk
        if pwi_params is None:
            pwi_params = self.pwi_params

        # set up the calculator
        calc = espresso.Espresso(pseudopotentials=self.pseudopotentials, tstress=False, tprnfor=False, kpts=nk,
                                 input_data=pwi_params)

        # get the mpi arguments -- these are all the arguments given to the program.
        mpi_args = ' '.join(sys.argv[1:])

        if mpi_args != '':
            calc.command = 'mpirun ' + mpi_args + ' pw.x -in PREFIX.pwi > PREFIX.pwo'

        # attach the calculator to the atoms
        atoms.set_calculator(calc)

        # calculate energy with DFT, and return energy
        try:
            energy = atoms.get_total_energy()
        except calculator.CalculationFailed:
            energy = None

        return energy

    def do_sweep(self, sweep_param_namespace, sweep_parameter, nk, parameters, atoms=None):
        if atoms is None:
            atoms = self.atoms

        successful_parameters = []
        e = []

        input_data = copy.copy(self.pwi_params)

        # for each parameter value...
        for parameter_value in parameters:
            # ... set up the input data
            try:
                # try to just write the new value into the dictionary - if it doesn't exist, add it
                input_data[sweep_param_namespace][sweep_parameter] = parameter_value
            except KeyError:
                # add the key and the value if not already in it
                input_data[sweep_param_namespace] = {sweep_parameter: parameter_value}

            # calculate the energy and add to the corresponding array
            energy = self.get_energy(atoms=atoms, nk=nk, pwi_params=input_data)

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

    def setup_dftmu(self, muon_positions, supercell: list=None, run_time=12*60*60, nodes=2,
                    devel=False, arc_priority_mode=False, prefix='PW', output_directory='dftmu', show_positions=True):
        """
        setup_dftmu: set up a DFT+mu calculation, all ready to sumbit on SLURM
        :param muon_positions: list of muon *scaled* positions (in unit cell) -- can generate using the utilities in
                               pointgrid.py, OR can give an integer, and will instead try to generate this number of
                               random positions (using pointgrid.py)
        :param run_time: run time in seconds
        :param nodes: number of nodes to run on
        :param supercell: list [nx,ny,nz] of how many unit cells in x y z to form the supercell. Default is [2,2,2]
        :param devel: make the slurm files to run in devel mode. Default is False
        :param arc_priority_mode: make the slurm files run in priority mode on ARC. Default is False (normal priority)
        :param prefix: prefix of the quantum espresso files. Suggested is [compoundname].relax.mu, but can be whatever
                       you want.
        :param output_directory: folder to store the files in.
        :param show_positions: if True, shows the muon positions
        :return: list of the names of files created
        """

        # if muon_positions is given as an integer, just make this many positions at random
        if isinstance(muon_positions, int):
            muon_positions = point_grid(input_structure=self.atoms,
                                        xgrid=muon_positions, ygrid=muon_positions, zgrid=muon_positions,
                                        random=muon_positions, rand_factor=muon_positions*2)

        # if supercell is not specified, do 2x2x2
        if supercell is None:
            supercell = [2, 2, 2]

        # define write slurm function
        def write_slurm(pw_file_name: str, directory=output_directory, runtime_s: int = run_time, n_nodes: int = nodes,
                        devel=devel, priority=arc_priority_mode, run_id=''):
            """
            Writes SLURM files to submit to the cluster
            :param pw_file_name: location of the pw.x file to run (if it doesn't have a .pwi extension it gets added)
            :param directory: directory on the cluster to save (and run) everything on
            :param runtime_s: runtime of the SLURM job, in seconds
            :param n_nodes: number of cores to run the script on
            :param devel: run on the development cores (for testing). If True it forces the runtime to be 600s.
            :param priority: run with priority quality of service on ARC
            :return: location of the SLURM file written
            """

            # sort out the file names - if it doesnt have the pwi extension then add it
            if pw_file_name[-4:] != '.pwi':
                pw_in_file_name = pw_file_name + '.pwi'
                pw_out_file_name = pw_file_name + '.pwo'
                slurm_file_name = pw_in_file_name + '_slurm.run'
            else:
                pw_in_file_name = pw_file_name
                pw_out_file_name = pw_in_file_name[:-4] + '.pwo'
                slurm_file_name = pw_in_file_name[:-4] + '_slurm.run'

            # if directory doesn't have a trailing /, add it
            if directory[-1] != '/':
                directory = directory + '/'

            # open the SLRUM file
            slurm_file = open(directory + slurm_file_name, 'w')

            # write the SLURM file
            slurm_file.write('#!/bin/bash\n\n')  # shebang
            slurm_file.write(
                '#======================================================================\n')  # start of parameters
            if priority:
                slurm_file.write('#SBATCH --qos=priority\n')
            # run on devel or compute partition?
            if devel:
                slurm_file.write('#SBATCH --partition=devel\n')
            else:
                slurm_file.write('#SBATCH --partition=compute\n')
            slurm_file.write('#SBATCH --constraint=haswell\n')
            slurm_file.write('#SBATCH --job-name=' + str(run_id) + pw_in_file_name[:-4] + '\n')  # job name
            if devel and runtime_s > 600:
                # if the runtime is beyond the devel partition limit, cap it
                slurm_file.write('#SBATCH --time=00:10:00\n')
            else:
                # otherwise just do as its told
                slurm_file.write('#SBATCH --time={:02d}'.format(int(np.floor(runtime_s/(60*60))))
                                 + time.strftime(":%M:%S", time.gmtime(runtime_s)) + '\n')
            slurm_file.write('#SBATCH --nodes=' + str(n_nodes) + '\n')  # number of nodes
            slurm_file.write('#SBATCH --ntasks-per-node=16\n')  # number of cores per node
            slurm_file.write(
                '#======================================================================\n\n')  # end of parameters

            # do modules
            slurm_file.write('module purge\nmodule load intel-compilers/2016\nmodule load intel-mkl/2016\n')
            slurm_file.write('module load mvapich2/2.1.0__intel-2016\nmodule load espresso/6.0.0\n\n')

            # init slurm environment variables
            slurm_file.write('. enable_arcus-b_mpi.sh\n')

            # move to the data directory
            slurm_file.write('cd $DATA/' + directory + '\n')
            # write out the pw.x command
            slurm_file.write('mpirun $MPI_HOSTS pw.x -i ' + pw_in_file_name + ' > ' + pw_out_file_name + '\n')

            # close the file
            slurm_file.close()

            # return the SLURM file location
            return slurm_file_name

        # convert the supercell into a matrix
        supercell = np.diag(supercell)

        supercell_atoms = build.make_supercell(self.atoms, supercell)

        if show_positions:
            all_muons_supercell = build.make_supercell(self.atoms, supercell)

        # give a hint about using the /scratch/ area for the output directory
        if self.pwi_params['control']['outdir'][0:8] != '/scratch':
            print('Hint: If running on ARC, it is better to use the /scratch area for the pw.x outdir.')

        # if 'relax' calcuation not specified, then warn
        try:
            if self.pwi_params['control']['calculation'] != 'relax':
                print('WARNING: Type of calculation is not set to relax. This is not normal DFT+mu!')
        except KeyError:
            print('You did not set the calculation type. No worries -- I\'m setting it to \'relax\' for you.')
            self.pwi_params['control'].update({'calculation': 'relax'})

        # write in the time limit
        self.pwi_params['control'].update({'max_seconds': run_time})
        try:
            self.pwi_params['control']['restart_mode'] = 'from_scratch'
        except KeyError:
            self.pwi_params['control'].update({'restart_mode': 'from_scratch'})

        # set up the calculator
        pw_calc = espresso.Espresso(pseudopotentials=self.pseudopotentials, tstress=False, tprnfor=True, kpts=self.nk,
                                    input_data=self.pwi_params)

        # for each muon site, append this to the crystal, and create the pw.x input file
        muon_inc = 0
        slurm_files = []
        pwi_files = []
        for muon_pos in muon_positions:
            # increment the counter for the prefix
            muon_inc += 1

            # convert the muon to angstrom (non-relative) coordinates
            # this is such a bodge! I'll sort it out when I have time
            muon_pos = np.asmatrix(self.atoms.get_cell()).transpose() * np.asmatrix(muon_pos).transpose()
            muon_pos_nparray = np.array(muon_pos)
            muon_pos = np.array([muon_pos_nparray[0][0], muon_pos_nparray[1][0], muon_pos_nparray[2][0]])

            # create the muon
            muon = atoms.Atom('H', muon_pos, mass=0.1134)

            # add the muon to the supercell
            supercell_atoms.append(muon)

            # sort out the prefix
            this_prefix = prefix + '.' + str(muon_inc)
            pw_calc.prefix = this_prefix
            pw_calc.label = output_directory + '/' + this_prefix  # label is basically the location of the file
            self.pwi_params['control']['prefix'] = this_prefix
            # write the pw.x input
            pw_calc.write_input(supercell_atoms)
            # write the SLURM input file
            slurm_files.append(write_slurm(pw_file_name=this_prefix, run_id=str(muon_inc)))
            pwi_files.append(this_prefix + '.pwi')

            # remove the muon from the supercell for the next iteration
            del supercell_atoms[-1]

            # add the muon onto the visual cell, so that it can be viewed
            if show_positions:
                all_muons_supercell.append(muon)

        # write the runall bash script to run on the cluster

        # if directory doesn't have a trailing /, add it
        if output_directory[-1] != '/':
            save_directory = output_directory + '/'

        script_location = output_directory + '/runall.sh'

        # create the runall file
        with open(script_location, 'w') as script_file:
            # do the shebang at the beginning
            script_file.write('#!/bin/bash\n')

            # for each slurm file...
            for slurm_file_name in slurm_files:
                # write the submission command
                script_file.write('sbatch ' + slurm_file_name + '\n')

        if show_positions:
            view(all_muons_supercell)

        return pwi_files

    def calculate_band_structure(self, k_points: list, n_inter_k_points=10, n_bands=None, cleanup=False):
        """
        Calculate the band structure using quantum espresso
        :param k_points: list of [[k point label, np.ndarray of k-point in conventional basis],...]
        :param n_inter_k_points: number of points between each of the pairs of k-points in k_points.
        :param n_bands: numbe of bands to calculate -- nbnd in pw.x
        :param cleanup: if True, will clean up all of the input/output files apart from the .dat and the .gle
        :return: True if OK, False if not
        """

        # set the output directory to 'out' if not set
        self.pwi_params = set_pwi_param(pwi_params=self.pwi_params, namespace='control', parameter_name='outdir',
                                        parameter_value='out')
        # set the calculation type to 'scf'
        self.pwi_params = set_pwi_param(pwi_params=self.pwi_params, namespace='control', parameter_name='calculation',
                                        parameter_value='scf', overwrite=True)
        # set the prefix to the compound name
        self.pwi_params = set_pwi_param(pwi_params=self.pwi_params, namespace='control', parameter_name='prefix',
                                        parameter_value=str(self.atoms.symbols))

        # first things first -- calculate the energy
        energy = self.get_energy()

        print('Calculated SCF energy is ' + str(energy) + ' eV')

        # do the k path
        k_path = []
        k_path_str = ""
        for i_k_point in range(0, len(k_points) - 1):
            # find the difference between this vector and the next one
            diff_k = k_points[i_k_point + 1][1] - k_points[i_k_point][1]
            # divide this difference by the number of points per line
            dk = diff_k / n_inter_k_points
            # and consecutively add this to k_points[i] and store in array and write to file
            k_j = k_points[i_k_point][1]
            for j in range(0, n_inter_k_points):
                k_path.append(k_j)
                k_path_str += "{:1.5f} {:1.5f} {:1.5f} 0 \n".format(k_j[0], k_j[1], k_j[2])
                k_j = k_j + dk
        k_path_str = str(len(k_path)) + '\n' + k_path_str

        # now we have the energy, we need to find the .pwi file -- which should be called espresso.pwi.
        with open('espresso.pwi', 'r') as pwi_file:
            # look for the lines that start with 'K_POINTS', 'calculation', '&SYSTEM',
            pwi_file_lines = pwi_file.readlines()

        calculation_type_line, k_points_line, system_line = None, None, None
        for line_id, line in enumerate(pwi_file_lines):
            if line.startswith('&SYSTEM'):
                system_line = line_id
            if line.lstrip().startswith('calculation'):
                calculation_type_line = line_id
            if line.startswith('K_POINTS'):
                k_points_line = line_id
                # k points will be the last thing -- so do this last
                break
        # now replace these lines with the new values
        # start with calculation:
        pwi_file_lines[calculation_type_line] = '   calculation      = \'bands\'\n'
        # now add in nbnd:
        pwi_file_lines.insert(system_line+1, '   nbnd     = ' + str(n_bands) + '\n')
        # inserting a line above means the K_POINTS line is now 1 more than it was before
        k_points_line += 1
        # now do the k points: get rid of the original two lines, replace with the k point path
        pwi_file_lines[k_points_line] = 'K_POINTS crystal\n'
        pwi_file_lines[k_points_line+1] = k_path_str

        with open('bands.pwi', 'w') as pwi_file:
            pwi_file.truncate()
            pwi_file.writelines(pwi_file_lines)

        # run the command for the band structure
        mpi_args = ' '.join(sys.argv[1:])
        if mpi_args != '':
            pw_bands_command = 'mpirun ' + mpi_args + ' pw.x -in bands.pwi > bands.pwo'
        else:
            pw_bands_command = "pw.x -in bands.pwi > bands.pwo"

        try:
            subprocess.run(pw_bands_command, shell=True, check=True)
        except subprocess.CalledProcessError:
            print("Sorry! The pw.x band structure calculation failed ☹️. Check the file bands.pwi and bands.pwo for "
                  "information")
            return False

        # looks like the pw.x bands calculation worked! Now do the bands.x calculation
        with open('bands.bandi', "w", encoding='UTF-8') as bands_file:
            bands_file.write("&BANDS\n" +
                             "   outdir   =  \'" + str(self.pwi_params['control']['outdir']) + "\',\n" +
                             "   prefix   =  \'" + str(self.pwi_params['control']['prefix']) + "\',\n" +
                             "   filband  =  \'" + str(self.pwi_params['control']['prefix']) + ".dat\',\n" +
                             "   lsym     =  .true,\n/\n")

        # now run bands.x
        if mpi_args != '':
            bands_command = 'mpirun ' + mpi_args + ' bands.x -i bands.bandi'
        else:
            bands_command = "bands.x -i bands.bandi"
        # we have to do it this way becasue of an annoying segmentation error (quantum espresso's fault!)
        bands_process = subprocess.run(bands_command, shell=True, check=False, universal_newlines=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bands_output = bands_process.stdout.split('\n')

        # get the high symmetry x-points
        high_symmetry_x = []
        print('These are the high symmetry points found by quantum espresso. Check they\'re correct!')
        try:
            for bands_output_line in bands_output:
                if 'high-symmetry point' in bands_output_line:
                    print(str(bands_output_line) + ' <--> ' + str(k_points[len(high_symmetry_x)][0]) + '..?')
                    high_symmetry_x.append(bands_output_line.split()[-1])
        except IndexError:
            print('❄️ Symmetry error. Check the symmetry points you gave are actually points of high symmetry.\n'
                  '(Don\'t worry -- you don\'t have to run all this again! Use plot_bands in gle_utils to do this '
                  'manually, if you want to proceed with these symmetry points anyway. All the DFT stuff has been done.'
                  ')')
            return False

        # if after all that it didn't find any high symmetry points, then complain
        if len(high_symmetry_x) == 0:
            print('No high symmetry points were found by bands.x 🤯')
            print('Proceeding anyway...')

        # find the highest occupied energy
        highest_occupied_energy = 0
        with open('espresso.pwo', 'r') as scf_file:
            scf_lines = scf_file.readlines()
            for scf_line in scf_lines:
                if 'highest occupied' in scf_line:
                    print(scf_line)
                    highest_occupied_energy = float(scf_line.split()[-1])

        # now get the .dat.gnu file
        dat_gnu_file_location = str(self.pwi_params['control']['prefix']) + '.dat.gnu'

        high_symmetry_points = [[high_sym_x, k_points[i_high_sym][0]] for i_high_sym, high_sym_x
                                in enumerate(high_symmetry_x)]

        # now plot the bands! (this replaces glebands.py)
        plot_bands(dat_gnu_file_location=dat_gnu_file_location,
                   symm_x_labels=high_symmetry_points,
                   fermi_energy=highest_occupied_energy)

        if cleanup:
            import os
            os.remove('espresso.pwi')
            os.remove('espresso.pwo')
            os.remove('bands.pwo')
            os.remove('bands.pwi')
            os.remove('bands.bandi')
            os.remove(str(self.pwi_params['control']['prefix']) + '.dat.rap')
            os.remove(str(self.pwi_params['control']['prefix']) + '.dat.gnu')

        return True


def set_pwi_param(pwi_params: dict, namespace: str, parameter_name: str, parameter_value, overwrite=False):
    """
    set_pwi_param: goes through pwi_params, and sets the parameter if not there, or (if overwrite=True) overwrites it
    :param pwi_params: dictionary of pwi_params
    :param namespace: namespace of parameter in question
    :param parameter_name: name of parameter in question
    :param parameter_value: value to set the parameter to
    :param overwrite: if True, ignores the value already set and sets to parameter_value anyway.
    :return: pwi_params with the parameter set
    """

    # see if the namespace is there, if not, add it
    try:
        pwi_params[namespace]
    except KeyError:
        pwi_params.update({namespace: {}})

    # see if the parameter is in the namespace
    try:
        if pwi_params[namespace][parameter_name] != '':
            # the parameter is here! so only set if overwrite is True
            if overwrite:
                print('Overwriting parameter ' + parameter_name + ' to ' + str(parameter_value) + '\n')
                pwi_params[namespace][parameter_name] = parameter_value
    except KeyError:
        pwi_params[namespace].update({parameter_name: parameter_value})

    return pwi_params