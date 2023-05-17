"""
pw.py -- unit to do pw.x calculations using a nice, clean, simple python interface

Created by John Wilkinson, 8/10/2020
"""

from ase import io, build, atoms
from ase.io import espresso as ioespresso
from ase.visualize import view
from ase.build import bulk
from ase.calculators import espresso, calculator
import copy
import numpy as np
import sys
import time
from .pointgrid import point_grid
from .gle_utils import plot_xy, plot_scatter, plot_bands
from .hubbard import Hubbard
import subprocess
import os


import functools
print = functools.partial(print, flush=True)


class PW(object):
    """
    PW class -- contains everything you could ever want to do with Quantum Espresso's pw.x utility
    """

    def __init__(self, pwi_params: dict, atoms: bulk, pseudopotentials: dict, nk: tuple, qe_version=None,
                 hubbard: Hubbard = None):
        print('☕️ Welcome to musr_espresso -- DFT+mu tools for Quantum Espresso\n')
        
        # get the mpi arguments -- these are all the arguments given to the program.
        command_args = sys.argv[1:]
        self.parallel = False
        self.pw_command = 'pw.x'
        self.mpi_command = 'mpirun'
        if '--pw_command' in command_args:
            pw_path_index = command_args.index('--pw_command')
            self.pw_command = command_args[pw_path_index+1]
            del command_args[pw_path_index+1]
            del command_args[pw_path_index]
        else:
            self.pw_command = 'pw.x'
        print(self.pw_command)
        # try to get the pw.x version from the command, if it's not given in the arguments
        if qe_version is None:
            try:
                pw_out = subprocess.run(self.pw_command, stdin=subprocess.DEVNULL, capture_output=True)
                for i_line, potential_v_line in enumerate(str(pw_out.stdout).split('\\n')):
                    if i_line > 20:
                        raise IndexError
                    line_components = potential_v_line.split()
                    try:
                        if line_components[0] == 'Program' and line_components[1] == 'PWSCF':
                            qe_version = line_components[2][2:]
                            break
                    except IndexError:
                        continue
                print('Found pw.x version ' + str(qe_version) + '. Not what you were expecting? Change the '
                      '--pw_command\'s directory (if using ezq, this is set in <musr_espresso>/config/ezq.conf).')
            except (IndexError, subprocess.SubprocessError):
                print('Can\'t get pw.x version, sorry. Write it in yourself when you create an instance of the '
                      'pw class. I\'ll still try to run anyway...')
        self.qe_version = qe_version

        if '--mpi_command' in command_args:
            self.parallel = True
            mpi_path_index = command_args.index('--mpi_command')
            self.mpi_command = command_args[mpi_path_index+1]
            print(self.mpi_command)
            del command_args[mpi_path_index+1]
            del command_args[mpi_path_index]

        self.mpi_args = ' '.join(command_args)
        if self.mpi_args == '':
            print('😴 Using serial mode. To run in parallel, just add -n <n_processors> to the end of the python '
                  'command')
        else:
            print('😁 Using parallel mode. MPI options given are ' + self.mpi_args + '\n')
        self.pwi_params = pwi_params
        self.hubbard = hubbard
        self.atoms = atoms
        self.pseudopotentials = pseudopotentials
        self.pp_names = [item for item, _ in self.pseudopotentials.items()]
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
        :param small_sf: scale factor to apply to the cell for the small cell energy calculation. Can be a tuple of
                         (xscale, yscale, zscale) *or* a float of the scale factor
        :param plot: do plot of results
        :return: list with entries [nk_id][param][dE]
        """

        # if small_sf is just one number, turn it in to a tuple of (x, y, z)
        if isinstance(small_sf, float):
            small_sf = (small_sf, small_sf, small_sf)

        # if nk is a list, then loop up it
        if isinstance(nk, list):
            result = []
            for i_nk, this_nk in enumerate(nk):
                # get one of the list elements, and re-call this function with *one* of the k-points, not plotting
                result.append(self.sweep_parameter(pwi_namespace, pwi_parameter, pwi_param_values, this_nk, plot=False))
            if plot:
                all_parameters = [res[0] for res in result]
                dE = [res[1] for res in result]
                plot_xy(all_parameters, dE, pwi_parameter + ' and scf energy, ' + str(self.atoms.symbols),
                        legend=["nk=" + str(i) for i in nk],
                        ytitle="$E(a)-E(({:.3f}, {:.3f}, {:.3f}).a)$ (Ry)".format(*small_sf),
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
                    xtitle=pwi_parameter,
                    ytitle="$E(a)-E(({:.3f}, {:.3f}, {:.3f}).a)$ (Ry)".format(*small_sf))

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

        # make the pw.x command
        if self.mpi_args != '' or self.parallel==True:
            calc.command = self.mpi_command + ' ' + self.mpi_args + ' '
        else:
            calc.command = ''

        calc.command += self.pw_command + " -in " + calc.template.inputname + " > " + calc.template.outputname
        print(calc.command)

        # attach the calculator to the atoms
        atoms.set_calculator(calc)

        energy = None
        # calculate energy with DFT, and return energy
        try:
            # i knowwww, we could just do atoms.get_total_energy(). BUT this doesn't work if you have multiple
            # hubbard U parameters (try it!) so we have to do everything manually instead...
            self.write_pw_input(atoms=atoms, nk=nk, pwi_params=pwi_params)
            # the below lines are the bits copied from ASE which run the calculation (I had to copy this
            # because otherwise the file will be rewritten)
            try:
                proc = subprocess.Popen(calc.command, shell=True, cwd=calc.directory)
            except OSError as err:
                # Actually this may never happen with shell=True, since
                # probably the shell launches successfully.  But we soon want
                # to allow calling the subprocess directly, and then this
                # distinction (failed to launch vs failed to run) is useful.
                msg = 'Failed to execute "{}"'.format(calc.command)
                raise EnvironmentError(msg) from err
            errorcode = proc.wait()
            if errorcode:
                raise calculator.CalculationFailed('Calculation failed with code.' + str(errorcode))
            calc.template.read_results()
            self.check_scf_accuracy()
            energy = calc.results['energy']
        except calculator.CalculationFailed:
            print('Calculation failed.')
            energy = None
        except AssertionError:
            # AssertionErrors in ASE are often due to it reading the file wrong  -- so see if there is a total energy
            # line, and just get that.
            print('ASE can\'t find the energy, so finding it manually...')
            with open('espresso.pwo', 'r') as out_file:
                for line in out_file:
                    split_line = line.split()
                    if len(split_line) <= 2:
                        continue
                    if split_line[0] == '!':
                        # this is the total energy line! so get it (and convert to Ry)
                        energy = float(split_line[-2]) * 13.60566
                        break
                if energy == None:
                    print('🤒 Oh no! I can\'t even get it from the file. Something\'s gone very wrong')
        return energy

    def check_scf_accuracy(self, pw_output_file='espresso.pwo'):
        """
        Checks the 'estimated scf accuracy' in the pw.x decreases on each iteration. If it seems to
        decrease and then suddenly increase again, its a sign of the lowest unoccupied and highest occupied
        states exchanging. According to the documentation, adding in a broadenign should fix this!
        :param pw_output_file: file to check for good convergence
        :return: True if everything is OK, False (with explanation) if there are issues.
        """
        with open(pw_output_file, 'r') as f:
            accuracies = []
            for line in f.readlines():
                if 'estimated scf accuracy' in line:
                    accuracies.append(float(line.split()[-2]))
                if 'End of self-consistent calculation' in line:
                    break

        accuracy_rel_diff = [(accuracies[i-1] - accuracies[i])/accuracies[i-1] for i in range(1, len(accuracies))]

        if any(accuracy_rel_diff) < 0:
            print('Warning -- SCF error decreases and then suddenly increases. Here are the accuracies (in Ry), '
                  'decide for yourself if you trust the results...')
            if any(accuracy_rel_diff[2:]) < -2:
                print('(there is a factor of two difference, so I wouldnt trust it if I were you... try adding '
                      'in some smearing.)')
            print(accuracies)

            return False
        return True

    def write_pw_input(self, atoms: bulk = None, nk = None, pwi_params = None, filename='espresso.pwi'):
        """
        Write the pw.x input file (mostly uses ASE's stuff, but also lets you do the element names instead
        of ids for element-dependent properties (e.g Hubbard U), which can crash ASE
        :param atoms: atoms to calculate on
        :param nk: k-grid to calculate on
        :param pwi_params:  pw.x input parameters
        :param filename: filename to write to
        :return:
        """

        # take out all the bits in the input parameters with brackets and put into a separate dict
        problem_keys = {}
        for namespace in self.pwi_params:
            current_namespace = {}
            for key in list(self.pwi_params[namespace]):
                if '(' in key:
                    current_namespace.update({key: self.pwi_params[namespace][key]})
                    del self.pwi_params[namespace][key]
            problem_keys.update({namespace: current_namespace})

        try:
            # try making the folder, but don't worry if it can't...
            os.makedirs(os.path.dirname(filename), exist_ok=True)
        except FileNotFoundError:
            pass
        with open(filename, 'w') as f:
            ioespresso.write_espresso_in(f, atoms, pwi_params, self.pseudopotentials,
                                         kpts=nk)

        # now put the problem keys into the pwi file:

        # read in the pwi file lines
        with open(filename, 'r') as f:
            lines = f.readlines()

        # find the pseudopotentials
        pp_i = lines.index('ATOMIC_SPECIES\n') + 1
        while lines[pp_i] != '\n':
            this_line_split = lines[pp_i].split()
            if not this_line_split[0] in self.pp_names:
                self.pp_names.append(this_line_split[0])
            pp_i += 1

        # for each problem_key namespace...
        for namespace in problem_keys:
            # ... find the line which is the beginning of the namespace
            namespace_line = '&' + namespace.upper() + '\n'
            namespace_line_id = lines.index(namespace_line)
            # now add in the problem bits
            for key in list(problem_keys[namespace]):
                species_identifier = key[key.find("(")+1:key.find(")")]
                if species_identifier in self.pp_names:
                    # the user used the atomic symbol instead of a numeric ID -- so see how many of these
                    # are in the pw.x file (there might be more than 1 for e.g AFM structures)
                    pp_index = self.pp_names.index(species_identifier)
                    while pp_index < len(self.pp_names):
                        if species_identifier in self.pp_names[pp_index]:
                            # we have a key <--> species match!
                            this_key = key.split('(')[0] + '(' + str(pp_index + 1) + ')'
                            lines.insert(namespace_line_id + 1, this_key + ' = ' + str(problem_keys[namespace][key]) + '\n')
                        pp_index += 1
                else:
                    # there doesn't seem to be a relation between the species identifier and the element names,
                    # so just put it in the file and don't faff
                    lines.insert(namespace_line_id+1, key + ' = ' + str(problem_keys[namespace][key]) + '\n')

        # now write it out again
        with open(filename, 'w') as f:
            file_contents = "".join(lines)
            f.write(file_contents)

            # and do the hubbard parameters if qe>=7
            if float(self.qe_version) > 7 and self.hubbard is not None:
                f.write(self.hubbard.qe_7_hubbard(pp_names=self.pp_names))

        # replace the problem keys in the dict
        for namespace in problem_keys:
            self.pwi_params[namespace].update(problem_keys[namespace])


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
                    devel=False, arc_priority_mode=False, prefix='PW', output_directory='dftmu',
                    cluster_directory="~", show_positions=True):
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
        :param cluster_directory: directory where the files will be on the cluster.
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
        def write_slurm(pw_file_name: str, cluster_directory=cluster_directory, directory=output_directory,
                        runtime_s: int = run_time, n_nodes: int = nodes, devel=devel, priority=arc_priority_mode,
                        run_id=''):
            """
            Writes SLURM files to submit to the cluster
            :param pw_file_name: location of the pw.x file to run (if it doesn't have a .pwi extension it gets added)
            :param cluster_directory: directory on the cluster to save (and run) everything on
            :param directory: directory to store the slurm files locally
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
            slurm_file.write('cd ' + cluster_directory + '\n')
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
            # pw_calc.write_input(supercell_atoms)
            self.write_pw_input(atoms=supercell_atoms, nk=self.nk, pwi_params=self.pwi_params,
                                filename=pw_calc.label + '.pwi')
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

        if self.mpi_args != '' or self.parallel==True:
            pw_bands_command = self.mpi_command + ' ' + self.mpi_args + ' ' + self.pw_command + ' -in bands.pwi > bands.pwo'
        else:
            pw_bands_command = self.pw_command + " -in bands.pwi > bands.pwo"

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

        bands_command = self.pw_command[:-4] + 'bands.x'
        # now run bands.x
        if self.mpi_args != '' or self.parallel==True:
            bands_command = self.mpi_command + ' ' + self.mpi_args + ' ' + bands_command + ' -i bands.bandi'
        else:
            bands_command = bands_command +  " -i bands.bandi"
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
