#!/usr/local/bin/python3.7
import argparse
import pandas as pd
from pathlib import Path
from ase import io, spacegroup


def main():

    # get the arguments
    parser = argparse.ArgumentParser(description='pw.x output file processor to find the muon site')
    parser.add_argument('-sg', '--spacegroup', nargs=1, default='1')
    parser.add_argument('-em', '--emax', nargs=1, required=False, type=float,
                        help="Maximum energy difference (in Kelvin) from the smallest run to put in CSV file")
    parser.add_argument('-c', '--cluster', action='store_true', help='Do only one muon from each cluster')
    parser.add_argument('-uc', '--unitcellcif', nargs=1, required=True, type=str,
                        help="CIF file of the original structure to put the final muon sites into")
    parser.add_argument('csv_file', type=str, help='CSV file output of pwprocess.py which has the details of '
                                                              'all the pw.x results')

    arguments = parser.parse_args()

    sg = arguments.spacegroup
    if isinstance(sg, list):
        sg = sg[0]
    if sg.isdigit():
        sg = int(sg)

    unit_cell_cif = arguments.unitcellcif
    if isinstance(unit_cell_cif, list):
        unit_cell_cif = unit_cell_cif[0]

    emax = arguments.emax
    if isinstance(emax, list):
        emax = float(emax[0])

    docluster = arguments.cluster
    pwprocess_output = arguments.csv_file

    newcif_location = Path(pwprocess_output).stem + '.cif'

    muon_labels, muon_positions = get_unit_cell_muon_positions(pwprocess_output=pwprocess_output, docluster=docluster,
                                                               emax=emax)

    make_csv(original_cif=unit_cell_cif, new_file_location=newcif_location,  muon_sites=muon_positions,
             muon_labels=muon_labels, sg=sg)

    return 0


def get_unit_cell_muon_positions(pwprocess_output: str, docluster: bool = False, emax: float = None) -> (list, list):
    """
    Get the unit cell muon positions from pwprocess.py's CSV output
    :param pwprocess_output: location of the output file form pwprocess.py
    :param docluster: if True, only do one muon from each dist cluster
    :param emax: maximum difference in energy between each structure and that of minimum energy to get the muon position
                 of (Kelvin).
    :return: (list of muon site labels, list of [x,y,z] fractional muon positions)
    """

    csv_data = pd.read_csv(pwprocess_output, error_bad_lines=False, warn_bad_lines=False)
    csv_data = csv_data.dropna()

    # if no clustering, find out what position the id number of the files is
    filename_label_pos = 0
    if len(csv_data) > 1 and not docluster:
        filename0_split = csv_data['file_name'][0].split('.')
        filename1_split = csv_data['file_name'][1].split('.')
        for i_element in range(0, len(filename0_split)):
            try:
                if filename0_split[i_element] != filename1_split[i_element]:
                    filename_label_pos = i_element
            except IndexError:
                filename_label_pos = 0

    muon_labels = []
    muon_sites = []

    current_cluster = 0
    for i_muonsite in range(0, len(csv_data)):
        this_label = ''
        # check the energy is within range
        if emax is not None:
            if float(csv_data['dE (K)'][i_muonsite]) > emax:
                break
        # check this is a new cluster (if applicable)
        if docluster:
            if csv_data['Dist cluster'][i_muonsite] == current_cluster:
                continue
            else:
                current_cluster = csv_data['Dist cluster'][i_muonsite]
                this_label = current_cluster
        else:
            # if ignoring clustering, get the label from the file names
            this_label = csv_data['file_name'][i_muonsite].split('.')[filename_label_pos]
        muon_labels.append(this_label)
        muon_unit_x = csv_data['Mu_unit_xf'][i_muonsite]
        muon_unit_y = csv_data['Mu_unit_yf'][i_muonsite]
        muon_unit_z = csv_data['Mu_unit_zf'][i_muonsite]
        muon_sites.append([muon_unit_x, muon_unit_y, muon_unit_z])

    return muon_labels, muon_sites


def make_csv(original_cif: str, new_file_location: str, muon_sites: list, muon_labels: list, sg=1) -> str:
    """
    make a CSV file of the
    :param original_cif: location of the CIF file of the original structure
    :param new_file_location: location of the CIF file to be created
    :param muon_sites: list of FRACTIONAL coordinates of the muon positions
    :param muon_labels: labels of the muon positions to use in the CIF file
    :param sg: spacegroup of the structure (uses this to get the other muon positions)
    :return: str of the file written to
    """

    # read the cif file
    unit_cell_atoms = io.read(original_cif)

    # make a new one
    io.write(new_file_location, unit_cell_atoms)

    # re-open the file
    new_cif = open(new_file_location, 'a+')
    print(muon_sites)
    for i_muon in range(0, len(muon_sites)):
        # find all the positions of this muon by symmetry

        muon_symmetry_equivalent = spacegroup.crystal(symbols=['X'], basis=[muon_sites[i_muon]], spacegroup=sg,
                                            cell=unit_cell_atoms.get_cell())
        muon_symmetry_equivalent = muon_symmetry_equivalent.get_scaled_positions()
        for muon_position in muon_symmetry_equivalent:
            # muon_position = muon_sym.get_scaled_positions()
            new_cif.write('  H' + str(muon_labels[i_muon]) + '\t1.0000 ' + str(muon_position[0]) + ' ' +
                          str(muon_position[1]) + ' ' + str(muon_position[2]) + ' Biso  1.0000  H\n')

    new_cif.close()


if __name__ == '__main__':
    main()
