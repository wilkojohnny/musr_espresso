#!/opt/anaconda3/envs/qe_env/bin/python3

# perturbed_distance_plotter.py - gets a pw.x output file and plots the perturbed distances of each of the
# atoms wrt the muon position - sort of like Fig 2 Moeller et al PRB 87, 121108(R) 2013
# John Wilkinson 7/8/2019

from ase import io  # for reading quantum espresso output files
import argparse  # for command line option parsing
from musr_espresso import gle_utils  # for plotting
import os
# import glob  # for listing out files

def main():

    # get the .pwo file from the command arguments
    parser = argparse.ArgumentParser(description='Plot a graph of perturbed distances of atoms from pw.x output')
    parser.add_argument("pwo_file", help="pw.x output file")
    parser.add_argument("pwi_file", help="pw.x input file (only needed if input atoms not in pwo file)",
                        default=None)
    parser.add_argument("--nogle", action="store_true", help="Don't make gle file for plot")
    arguments = parser.parse_args()

    do_gle = not arguments.nogle

    # load the input file in
    pwi_file = arguments.pwi_file
    pwo_file = arguments.pwo_file

    output_prefix = os.path.basename(pwo_file)[:-4]

    dat_file_location = output_prefix + '.dist'
    gle_file_location = output_prefix + '.gle'

    # symbol used for the muon
    muon_symbol = 'H'

    # get the initial coordinates of the atoms
    if pwi_file is None:
        # if the user hasn't specified a pwi, file, get the initial positions from the pwo file
        initial_positions = io.read(pwo_file, index=0)
    else:
        initial_positions = io.read(pwi_file)

    # get the final coordinates of the atoms
    final_positions = io.read(pwo_file)

    # find the final position of the muon
    muon = None
    final_muon_id = None
    for atom in final_positions:
        if atom.symbol == muon_symbol:
            final_muon_id = atom.index
            muon = atom
            break

    # check that a muon has been found
    assert muon is not None

    # replace the muon in the initial_positions with the final location
    initial_muon_id = None
    for atom in initial_positions:
        if atom.symbol == muon_symbol:
            initial_muon_id = atom.index
            initial_positions[initial_muon_id].position = muon.position

    # calculate the minimum difference between each of the atoms and the muon:
    #    (minimum in the sense that +- lattice translation vectors until we get something small)
    distances = initial_positions.get_distances(initial_muon_id, range(0, initial_positions.get_number_of_atoms()),
                                                mic=True)
    distances = [distances, final_positions.get_distances(final_muon_id,
                                                          range(0, final_positions.get_number_of_atoms()), mic=True)
                 - distances]

    # swap around the array indices
    final_distances = []
    for i in range(0, len(distances[0])):
        final_distances.append([final_positions[i].symbol, distances[0][i], distances[1][i]])
    final_distances.sort(key=lambda x: (x[0], x[1]))

    # write the data to a .dat file
    current_species = ''
    all_species = []
    dat_file = None
    dat_file_locations = []
    for atom in final_distances:
        # don't bother plotting the muon
        if atom[0] != muon_symbol:
            # if the current species is a new one, make a new file for it
            if current_species != atom[0]:
                current_species = atom[0]
                all_species.append(current_species)
                if dat_file is not None:
                    dat_file.close()
                dat_file_locations.append(dat_file_location + '.' + current_species + '.dat')
                dat_file = open(dat_file_locations[-1], 'w')
                dat_file.write('! unperturbed_distances_ang\tperturbed_distances_ang\n')
            if dat_file is not None:
                dat_file.write(str(atom[1]) + '\t' + str(atom[2]) + '\n')

    if do_gle:
        gle_utils.plot_scatter(dat_file_locations, title=pwo_file, legend=all_species, file_location=gle_file_location)

    return


if __name__ == '__main__':
    main()
