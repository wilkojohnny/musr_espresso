# qe_inputwriter.py - batch write input files for DFT calculations
# John Wilkinson 9/5/19

from ase.calculators.espresso import Espresso  # for quantum espresso calculations
from ase.visualize import view  # for visualising the crystal structure
from ase import io, build, atoms  # for reading in CIF files, to build supercells, and individual atoms
import numpy as np  # for matrices and numerical stuff
import pointgrid  # Franz's script to create muon starting locations (modified to allow premade structures)
import os  # for making directories
import time  # for converting time formats

# standard inputs which work generally
input_data = {'CONTROL': {
                            'calculation': 'relax',
                            'pseudo_dir': 'pp',
                            'outdir': 'out',
                            },
              'IONS': {
                        'ion_dynamics': 'bfgs'
                      },
              'SYSTEM': {},
              'ELECTRONS': {}
              }

muon_mass = 0.1134


def main():
    # define pseudopotential file names
    pseudopotentials = {'Na': 'Na.pbe-spnl-rrkjus_psl.1.0.0.UPF',
                        'P': 'P.pbe-nl-rrkjus_psl.1.0.0.UPF',
                        'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
                        'F': 'F.pbe-n-rrkjus_psl.0.1.UPF',
                        'H': 'H.pbe-rrkjus_psl.1.0.0.UPF'}

    nk = (5, 5, 3)  # number of k points
    prefix = 'Na2PO3F.relax.mu'  # prefix for pw.x
    runtime = 24*60*60 - 1  # maximum runtime of pw.x
    ciffile_location = '~/Documents/University/Na2PO3F/Na2PO3F.cif'  # location of .cif structure file
    output_directory = 'Na2PO3F'
    n_nodes = 5

    # define additional inputs
    input_data['CONTROL'].update({'prefix': prefix,
                                  'restart_mode': 'from_scratch',
                                  'max_seconds': runtime})
    input_data['SYSTEM'] = {'tot_charge': +1.0,
                            'ecutwfc': 85.,
                            'ecutrho': 900.,
                            'nspin': 1,
                            # 'lda_plus_u': True,
                            # 'Hubbard_U(1)': 3,
                            # 'occupations': 'smearing',
                            # 'smearing': 'gaussian',
                            # 'degauss': 0.01
                            }
    input_data['ELECTRONS'] = {'mixing_beta': 0.7}

    # make the output directory if it doesn't already exist
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # get the crystal cell
    crystal_cell = io.read(ciffile_location)

    # use pointgrid to find the muon positions
    # muon_positions = pointgrid.get_internuclei_positions(crystal_cell, ['F'], 4, spg=221, prec_grid=0.2)

    muon_positions = pointgrid.point_grid(crystal_cell, 5, 5, 5, random=15, rand_factor=10, set_spg=19)

    # construct the supercell
    crystal_supercell = build.make_supercell(crystal_cell, np.diag([1, 1, 1]))

    # set up the calculator
    pw_calc = Espresso(pseudopotentials=pseudopotentials, tstress=False, tprnfor=True, kpts=nk,
                       input_data=input_data)

    # for each muon site, append this to the crystal, and create the pw.x input file
    muon_inc = 0
    slurm_files = []
    for muon_pos in muon_positions:
        # increment the counter for the prefix
        muon_inc += 1

        # convert the muon to angstrom (non-relative) coordinates
        # this is such a bodge! I'll sort it out when I have time
        muon_pos = np.asmatrix(crystal_cell.get_cell()).transpose() * np.asmatrix(muon_pos).transpose()
        muon_pos_nparray = np.array(muon_pos)
        muon_pos = np.array([muon_pos_nparray[0][0], muon_pos_nparray[1][0], muon_pos_nparray[2][0]])

        # create the muon
        muon = atoms.Atom('H', muon_pos, mass=0.1134)

        # add the muon to the supercell
        crystal_supercell.append(muon)

        # sort out the prefix
        this_prefix = prefix + '.' + str(muon_inc)
        pw_calc.prefix = this_prefix
        pw_calc.label = output_directory + '/' + this_prefix  # label is basically the location of the file
        input_data['CONTROL']['prefix'] = this_prefix
        # write the pw.x input
        pw_calc.write_input(crystal_supercell)
        # write the SLURM input file
        slurm_files.append(write_slurm(this_prefix, output_directory, runtime_s=runtime, n_nodes=n_nodes, devel=False,
                                       run_id=str(muon_inc)))

        # remove the muon from the supercell for the next iteration
        del crystal_supercell[-1]

        # add the muon onto the single cell, so that it can be viewed
        crystal_cell.append(muon)

    view(crystal_cell)

    # write the runall bash script to run on the cluster
    write_runall_script(slurm_files, output_directory)

    return 1


def write_slurm(pw_file_name: str, directory='', runtime_s: int = 43200, n_nodes: int = 1, devel=False, run_id=''):
    """
    Writes SLURM files to submit to the cluster
    :pw_file_name: location of the pw.x file to run (if it doesn't have a .pwi extension it gets added)
    :directiry: directory on the cluster to save (and run) everything on
    :runtime_s: runtime of the SLURM job, in seconds
    :n_cores: number of cores to run the script on
    :devel: run on the development cores (for testing). If True it forces the runtime to be 600s.
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
    slurm_file.write('#======================================================================\n')  # start of parameters
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
        slurm_file.write('#SBATCH --time=' + time.strftime("%H:%M:%S", time.gmtime(runtime_s)) + '\n')
    slurm_file.write('#SBATCH --nodes=' + str(n_nodes) + '\n')  # number of nodes
    slurm_file.write('#SBATCH --ntasks-per-node=16\n')  # number of cores per node
    slurm_file.write('#======================================================================\n\n')  # end of parameters

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


def write_runall_script(slurm_file_names: list, save_directory=''):
    """
    Write a BASH script which submits all the generated SLURM files to the cluster
    :slurm_file_names: list of file names of the SLURM files to run
    :save_directory: directory to save the output file to
    :return: script_location (str)
    """

    # if directory doesn't have a trailing /, add it
    if save_directory[-1] != '/':
        save_directory = save_directory + '/'

    script_location = save_directory + '/runall.sh'

    # create the file
    script_file = open(script_location, 'w')

    # do the shebang at the beginning
    script_file.write('#!/bin/bash\n')

    # for each slurm file...
    for slurm_file_name in slurm_file_names:
        # write the submission command
        script_file.write('sbatch ' + slurm_file_name + '\n')

    # close the file
    script_file.close()

    # return the location of the script file
    return script_location


if __name__ == '__main__':
    main()
