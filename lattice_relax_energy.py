#!/usr/local/bin/python3.7
"""
lattice_relax_energy.py -- from a completed pw.x relaxed calculation, interpolate N points between
unrelaxed and relaxed lattice, finding the energy of all of them.
John Wilkinson 1/9/20
"""

import ase
import ase.io
import qe_inputwriter
import os
import argparse


def main():
    parser = argparse.ArgumentParser(description='make pw.x input files from unrelaxed structure to relaxed structure,'
                                                 'to calculate the energy profile as the structure relaxes')
    parser.add_argument('-p', '--prefix', required=True, type=str)
    parser.add_argument('-N', '--nfiles', required=True, type=int,
                        help="Number of steps between relaxed and unrelaxed")
    parser.add_argument('-nn', '--nnodes', type=int, help='Number of nodes on the cluster to use')

    arguments = parser.parse_args()

    # number of steps
    N = arguments.nfiles - 1

    # get the pw.x input and output files
    file_prefix = arguments.prefix

    n_nodes = arguments.nnodes

    directory = str.split(file_prefix, '.')[0] + '_latticerelax'
    os.makedirs(directory, exist_ok=True)

    # copy all the preamble from the .pwi file into memory -- this will be the basis of the new file
    input_file = open(file_prefix + '.pwi', 'r')
    pwi_preamble = ''
    for pwi_line in input_file.readlines():
        if pwi_line.startswith('ATOMIC_POSITIONS'):
            break
        elif '\'relax\'' in pwi_line:
            # change the calculation type sneakily
            pwi_preamble += '   calculation = \'scf\'\n'
        else:
            pwi_preamble += pwi_line

    input_file.close()

    # find the old and new atomic positions
    initial_atoms = ase.io.read(file_prefix + '.pwi')
    final_atoms = ase.io.read(file_prefix + '.pwo', -1)
    muon = None

    # get the muon location
    for i_atom, atom in enumerate(final_atoms):
        if atom.symbol == 'H':
            initial_atoms.pop(i_atom)
            muon = final_atoms.pop(i_atom)

    slurm_files = []

    # for i=0 (unrelaxed lattice) to 10 (fully relaxed lattice), make a pw.x input file
    for i in range(0, N+1):
        this_file_name = directory + '/' + file_prefix + '.{:d}.pwi'.format(i)
        this_preamble = str.replace(pwi_preamble, file_prefix, file_prefix + '.{:d}'.format(i), 1)
        this_pwi_file = open(this_file_name, 'w+')
        # write in the preamble
        this_pwi_file.write(this_preamble)

        # write in the atoms
        this_pwi_file.write('\nATOMIC_POSITIONS angstrom\n')
        for i_atom, final_atom in enumerate(final_atoms):
            initial_position = initial_atoms[i_atom].position
            final_position = final_atom.position
            this_position = initial_position + (final_position - initial_position)*i/N
            this_pwi_file.write(final_atom.symbol +
                                ' {:12.8f} {:12.8f} {:12.8f}\n'.format(this_position[0], this_position[1], this_position[2]))
        this_pwi_file.write('H {:12.8f} {:12.8f} {:12.8f}\n'.format(muon.position[0], muon.position[1], muon.position[2]))
        this_pwi_file.close()
        slurm_files.append(qe_inputwriter.write_slurm(str.split(this_file_name[:-4], '/')[-1], directory, n_nodes=n_nodes,
                                                      devel=False, run_id=i))
    qe_inputwriter.write_runall_script(slurm_files, save_directory=directory)

    return 0


if __name__ == '__main__':
    main()
