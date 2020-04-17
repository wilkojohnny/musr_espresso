# -*- coding: utf-8 -*-
from __future__ import print_function  # python 2-3 compatible printing function
from ase import io, Atoms, units, spacegroup  # https://wiki.fysik.dtu.dk/ase/ase/atoms.html#list-of-all-methods
import numpy as np
import pandas as pd
import csv, glob, os, sys, operator, copy

# sys.path.insert(0, os.path.abspath('..')) # This to add the Soprano path to the PYTHONPATH so we can load it without installing it
sys.path.insert(0, os.path.abspath("/Users/johnny/Documents/University/DFT/qe_utils/Franz_DFT_scripts"))
from soprano.collection import AtomsCollection
from soprano.analyse.phylogen import Gene, PhylogenCluster
from utilities import utilities


def QErel(prefix, postfix='.out', set_spg=None, set_spgsetting=None, cif_flag=True, XMufile_flag=False, Esort_flag=True,
          set_ncluster=None, muon_multipl_prec=0.01, dist_cluster_prec=0.1, supercell=None, verbose='min'):
    """A function to analyse a set of Quantum Espresso output files for DFT+mu calculations Frang Lang 03/2018

    | Args:
    |   prefix (str): common prefix to all output files (has to include directory path if files are not in current folder)
    |   postfix (str): common postfix to all output files
    |   NB: all files of the form prefix*postfix will be used
    |   set_spg (str or int): spagegroup of unperturbed material. Can be in H-M notation or as number in Int. Tables. If not set: spglib will be used to determine it.
    |   set_spgsetting (int): setting of spacegroup, if an unconventional setting is needed. If not set it will revert to default value of 1.
    |   cif_flag (bool,default=True): whether to create a cif file with the final structure from the QE output
    |   XMufile_flag (bool,default=False): flag whether to create a copy of the QE output file with Mu->XMU and those new files will be used. If false, the original file will be changed only.
    |   Esort_flag (bool,default=True): whether the final results should be sorted in ascending energy of the final configuration.
    |   set_ncluster (int): number of clusters for the Phylogen clustering analysis. If not defined it will be set to the number of cluster found in the muon-muon distance clustering.
    |   muon_multipl_prec (float): Minimum 'distance' in scaled coordinates to determine multiplicity of each final muon site
    |   dist_cluster_prec (float): distance between final muon sites (in Angstrom) below which they are classed as belonging to the same cluster
    |   supercell (list of int): supercell used to calculate the relaxation.
    |   verbose (str: min,max,progbar, default=min): sets in how much detail the program reports on its progress. If no progress is wanted set it to any other string not specified in the list.
    """

    # sort out supercell input - if not given assume 1x1x1 unit cells
    if supercell is None:
        supercell = [1, 1, 1]

    filelist = glob.glob(
        prefix + '*' + postfix)  # find all files with names prefix*postfix (prefix has to include the directory if the files are not in the current folder)
    N_files = len(filelist)  # number of files found
    if verbose in ['min', 'max']:
        print(str(N_files), ' files found.')
    for filename in filelist:
        if XMufile_flag:
            utilities.file_strchange(filename, 'Mu', 'XMU',
                                     newfile_ext='XMU')  # change all occurrences of 'Mu' to 'XMU' so that ASE can read in the QE file (XMU is used to avoid turning XMu into XXMu etc.)
            if verbose in ['max']:
                print('Created a QE output file copy to: ', filename.replace(postfix, '.XMU' + postfix))
        else:
            utilities.file_strchange(filename, 'Mu', 'XMU')  # this one does NOT create a copy of the QE output file
    if XMufile_flag:  # add '.XMU' to the filename if XMufile_flag was set so that subsequently the copied QE files are used
        filelist = [x.replace(postfix, '.XMU' + postfix) for x in filelist]

    # create some arrays to fill later
    # list for column headers of the results summary
    results_headers = ['Mu_xi', 'Mu_yi', 'Mu_zi', 'Mu_xf', 'Mu_yf', 'Mu_zf', 'E (Ry)', 'dE (Ry)', 'dE (K)', 'NN ion',
                       'NN dist (A)', 'Multipl', 'Dist cluster', 'Bond order cluster', 'file_name']
    results_detailedheaders = ['Atom', 'x_i', 'y_i', 'z_i', 'x_f', 'y_f', 'z_f', 'Mu-dist.(A)', 'Displ.(A)', 'Cart_x_i',
                               'Cart_y_i', 'Cart_z_i', 'Cart_x_f', 'Cart_y_f', 'Cart_z_f',
                               'Change in dist. to final muon site (A)']
    results_summary = pd.DataFrame(index=range(0, N_files), columns=results_headers)  # array for the results summary
    resultsf_all = []  # array for all final configurations from relaxation calculations
    resultsi_all = []  # array for all initial configurations from relaxation calculations
    units2006 = units.create_units(
        '2006')  # raised the issue of espresso files being read with values for physical constants from 2006 with developers. Might be changed at some point soon.

    # look at each file in turn
    for iterator, filename in enumerate(filelist):
        if verbose in ['min', 'progbar']:  # write a progress bar
            sys.stdout.write(
                "\rFile: " + str(iterator + 1).rjust(len(str(N_files))) + '/' + str(N_files) + ' {0}'.format(
                    utilities.progbar(iterator + 1, N_files)))
        if verbose in ['max']:
            print('Reading file: ', filename)
        file_f = io.read(filename)  # read final (default) configuration
        file_i = io.read(filename, index=0)  # read initial configuration
        resultsf_all.append(file_f)  # add recently read final configuration to array holding all final configurations
        resultsi_all.append(
            file_i)  # add recently read initial configuration to array holding all initial configurations

        # write cif file of final configuration if cif_flag is set
        if cif_flag:
            file_f.write(filename.replace(postfix, '.cif'),
                         format='cif')  # this is the ASE write command for ASE atoms objects, NOT the default python write function
            if verbose in ['max']:
                print('Final configuration written to:', filename.replace(postfix, '.cif'))
        # after reading the first file, find number of atoms to create arrays to store the detailed information
        if iterator == 0:
            N_atoms = file_f.get_number_of_atoms()
            results_detailedpandas = np.array(np.zeros((N_files, N_atoms + 2, len(results_detailedheaders))),
                                              dtype=str)  # have two extra rows for: column headers and for summary row

        atoms_symbols = file_i.get_chemical_symbols()  # list of atom labels
        fractcoords_i = file_i.get_scaled_positions()  # initial fractional coordinates
        fractcoords_f = file_f.get_scaled_positions()  # final fractional coordinates
        cartcoords_i = file_i.get_positions()  # initial Cartesian coordinates
        cartcoords_f = file_f.get_positions()  # final Cartesian coordinates
        final_energy = file_f.get_total_energy() / units2006[
            'Rydberg']  # total energy of final scf calculation in Rydberg
        cartcoords_change = cartcoords_f - cartcoords_i  # vector change in Cartesian position of each atom (in Angstrom)
        cartcoords_dist = np.sqrt((cartcoords_change * cartcoords_change).sum(
            axis=1))  # distance (in Angstrom) from final to initial position of each atom. This appears to be much faster than using linal.norm, see: https://stackoverflow.com/questions/9171158/how-do-you-get-the-magnitude-of-a-vector-in-numpy

        results_summary.at[iterator, ['file_name']] = filename
        results_detailedpandas[
            iterator, range(2, N_atoms + 2), 0] = atoms_symbols  # list of chemical symbols for the ions
        results_summary.at[iterator, ['Mu_xi', 'Mu_yi', 'Mu_zi']] = np.round(fractcoords_i[-1],
                                                                             3)  # initial muon position. rounded to 3 decimals to remove floating point glitches
        results_detailedpandas[iterator, range(2, N_atoms + 2), 1:4] = np.array(np.round(fractcoords_i, 3),
                                                                                dtype=str)  # save all initial fractional positions
        results_summary.at[iterator, ['Mu_xf', 'Mu_yf', 'Mu_zf']] = fractcoords_f[-1]  # final fractional muon position
        results_detailedpandas[iterator, range(2, N_atoms + 2), 4:7] = np.array(fractcoords_f,
                                                                                dtype=str)  # save all final positions
        Mu_nndist = file_f.get_distances(-1, range(0, N_atoms),
                                         mic=True)  # calculate all (final) distances (in Angstrom) to the muon (accounting for equivalent atoms in adjacent cells)
        results_detailedpandas[iterator, range(2, N_atoms + 2), 7] = np.array(Mu_nndist, dtype=str)
        results_detailedpandas[iterator, range(2, N_atoms + 2), 8] = np.array(cartcoords_dist, dtype=str)
        results_detailedpandas[iterator, range(2, N_atoms + 2), 9:12] = np.array(cartcoords_i, dtype=str)
        results_detailedpandas[iterator, range(2, N_atoms + 2), 12:15] = np.array(cartcoords_f, dtype=str)

        results_summary.at[iterator, 'E (Ry)'] = final_energy  # total energy of final scf calculation in Rydberg

        NN_mindist_index, NN_mindist_value = min(enumerate(Mu_nndist[0:-1]), key=operator.itemgetter(
            1))  # find the smallest final distance to the muon and the index of that atom
        results_summary.at[iterator, 'NN ion'] = atoms_symbols[
            NN_mindist_index]  # chemical symbol of muon's nearest neighbour ion
        results_summary.at[iterator, 'NN dist (A)'] = NN_mindist_value  # nearest neighbour ion distance to muon

        Mu_nndistchange = copy.deepcopy(
            file_i)  # create an actual copy of the initial configuration. Necessary as otherwise the next line will change file_i!
        Mu_nndistchange.pop(i=-1)  # remove initial muon site from initial configuration
        Mu_nndistchange.append(file_f[-1])  # add final muon site to initial configuration
        Mu_nndistchange = Mu_nndistchange.get_distances(-1, range(0, N_atoms),
                                                        mic=True)  # distance between all initial atom positions and final muon position
        results_detailedpandas[iterator, range(2, N_atoms + 2), 15] = np.array(Mu_nndist - Mu_nndistchange,
                                                                               dtype=str)  # change between distances to final muon position

        # write first two rows for each run details: short summary and then column headers.
        results_detailedpandas[iterator, 1, range(len(results_detailedheaders))] = np.array(results_detailedheaders,
                                                                                            dtype=str)
        results_detailedpandas[iterator, 0, range(len(results_detailedheaders))] = np.array(
            [atoms_symbols[-1], fractcoords_i[-1, 0], fractcoords_i[-1, 1], fractcoords_i[-1, 2], fractcoords_f[-1, 0],
             fractcoords_f[-1, 1], fractcoords_f[-1, 2], NN_mindist_value, 'Final', 'Energy', final_energy, 'Ry',
             final_energy * units.Ry, 'eV', final_energy * units.Ry / units.kB, 'K'], dtype=str)

    results_summary.at[:, 'dE (Ry)'] = results_summary.loc[:, 'E (Ry)'] - np.amin(
        results_summary.loc[:, 'E (Ry)'])  # subtract minimum energy from all other energies
    results_summary.at[:, 'dE (K)'] = results_summary.loc[:, 'dE (Ry)'] * units.Ry / units.kB
    if verbose in ['min', 'max']:
        print('\nQE output read in and summarised.')

    # determining the spacegroup of the unperturbed system at hand through spglib
    original_latt = copy.deepcopy(
        file_i)  # load the initial configuration from the last file read in. NEEDS the deepcopy function as otherwise the .pop() function will change the original!
    original_latt.pop(
        i=-1)  # Remove and return atom at index i (default last), which here is the muon. (to do: generalise this)
    spg, num_syms = utilities.get_spg(original_latt, num_syms_flag=True, verbose=False)
    if set_spg is None:
        if verbose in ['min', 'max']:
            print('Spacegroup detected from file: ', str(spg))
        spg = spg.split()[0]  # H-M of spacegroup only. Element [1] would be number in international tables
    else:
        spg = set_spg
        if verbose in ['min', 'max']:
            print('Spacegroup set by user: ', str(spg))
    if set_spgsetting is None:
        spgsetting = 1
        if verbose in ['min', 'max']:
            print('Default setting: 1')
    else:
        spgsetting = set_spgsetting
        if verbose in ['min', 'max']:
            print('User setting: ', str(spgsetting))

    # analysis of final muon sites: calculate all the muon-muon distances in an N_files x N_files array
    muon_distances = np.zeros((N_files, N_files))  # create an empty array for the muon-muon distances to go into
    muon_multiplicities = np.zeros(N_files)  # empty list for multiplicy of each final muon position
    muon_symbols = [
                       'X'] * N_files  # list of ['X','X',...] with number of elements equal to number of DFT+mu relaxation calculations
    muon_finalsites = results_summary.loc[:, ['Mu_xf', 'Mu_yf',
                                              'Mu_zf']].values  # list of fractional coordinates of final relaxed muon sites. as_matrix() is needed for ASE to recognize it as an array.
    # <JW>: change the muon sits so that they are all in the same unit cell of the supercell by doing (frac_coord % 1/supercellnes) * supercellness
    muon_finalsites[:, 0] = muon_finalsites[:, 0] % (1 / supercell[0]) * supercell[0]
    muon_finalsites[:, 1] = muon_finalsites[:, 1] % (1 / supercell[1]) * supercell[1]
    muon_finalsites[:, 2] = muon_finalsites[:, 2] % (1 / supercell[2]) * supercell[2]
    cell = np.diag(1/np.array(supercell))*original_latt.get_cell()
    # </JW> (also changed original_latt.getcell() to cell defined above)
    muon_finalatoms = Atoms(muon_symbols, muon_finalsites,
                            cell=cell)  # create Atoms object with the final relaxed muon sites. This will later be added to a lattice corresponding to one muon site in the original spacegroup
    muon_finalatoms.set_scaled_positions(
        muon_finalsites)  # set scaled positions again since ASE for some reason will not use scaled positions automatically
    for iterator in range(0, N_files):
        muon_latt = spacegroup.crystal(symbols=muon_symbols[iterator], basis=muon_finalsites[iterator], spacegroup=spg,
                                       setting=spgsetting, cell=cell, symprec=muon_multipl_prec)
        # this creates a lattice from the original (or passed) spacegroup and unit cell dimensions and with one of the final muon sites as the basis
        muon_multiplicities[
            iterator] = muon_latt.get_number_of_atoms()  # multiplicity of final muon site, up to symmetry precision set in previous line
        muon_latt.extend(
            muon_finalatoms)  # add all final muon sites to created lattice, in order to calculate the shortest distance between all symmetry equivalent positions of a particular muon site with all others (including itself)
        muon_distances[iterator] = np.amin(
            muon_latt.get_all_distances(mic=True)[-N_files:, 0:np.int_(muon_multiplicities[iterator])], axis=1)
        # get the distances from the symmetry equivalent positions of the currently selected final muon site to all other final muon sites (hence only select the last N_files rows and the first number of columns corresponding to the number of symmetry equivalent sites of currently selected final muon site), then select the minimum distance for each one.
        results_summary.at[iterator, 'Multipl'] = muon_multiplicities[
            iterator]  # add multiplicities to the results summary
    if verbose in ['min', 'max']:
        print('Distances between final muon sites calculated.')

    # clustering analysis purely based on muon-muon distances
    cluster_counter = 0  # counter for counting the number of clusters have been found
    dist_cluster = np.zeros(N_files)  # empty list to store the label of the cluster that each muon site belongs to
    for irow in range(0, N_files):  # cycle through the final muon sites
        if dist_cluster[irow] == 0:  # if the currently considered muon site has not been assigned a cluster then ...
            cluster_counter = cluster_counter + 1  # ... increment the cluster label by 1 ...
            for icolumn in range(0, N_files):  # ... cycle through all the columns (in the given row) and ...
                if muon_distances[
                    irow, icolumn] <= dist_cluster_prec:  # ... if the muon-muon distance is less than the precision set ...
                    dist_cluster[icolumn] = cluster_counter  # ... assign the cluster label to that muon.
        results_summary.at[irow, 'Dist cluster'] = dist_cluster[irow]  # add distance clustering to results summary
    if verbose in ['min', 'max']:
        print('Distance clustering done. Found ', str(cluster_counter), ' clusters.')

    # clustering analysis using soprano
    # first decide what the number of clusters should be
    if set_ncluster is None:
        ncluster = cluster_counter  # number of clusters from the distance clustering
    else:
        ncluster = set_ncluster  # externally set number of clusters
    aColl = AtomsCollection(filelist, progress=False)  # load in all output files
    # gene_lnk=Gene(name='linkage_list', weight=1.0, params={})#genes based on all atom pair distances
    # gene_E=Gene(name='energy',weight=1.0,params={})#energy from the ASE get_potential_energy() function
    gene_bonds = Gene(name='bond_order_pars', weight=1.0, params={
        'cutoff_radius': 3.0})  # bond order parameters, which are projections of the local environment of an atom on spherical harmonics
    phClust = PhylogenCluster(aColl, [gene_bonds])  # create the cluster for the given choice of gene(s)
    clust_inds, clust_slices = phClust.get_kmeans_clusters(
        ncluster)  # Returns indices and slices representing the clusters
    results_summary.at[
        range(0, N_files), 'Bond order cluster'] = clust_inds  # add clustering to the results summary table
    if verbose in ['min', 'max']:
        print('Bond order parameter clustering done with ', str(ncluster), ' clusters.')

    # sorting results by final energy
    index_Esorted = [[x, results_summary.loc[x, 'E (Ry)']] for x in
                     range(0, N_files)]  ##array with index and energy for the configuration with that index
    if Esort_flag:
        sort_column = 1  # sort by energy
    else:
        sort_column = 0  # sort by index (i.e. unsorted)
    index_Esorted = sorted(index_Esorted, key=lambda x: x[sort_column])
    # rejuggling the muon_distances array based on the sorted muon sites
    muon_distances_sorted = muon_distances
    for i in range(0, N_files):
        for j in range(0, N_files):
            muon_distances_sorted[i, j] = muon_distances[index_Esorted[i][0], index_Esorted[j][0]]
    muon_distances = muon_distances_sorted
    if verbose in ['min', 'max']:
        print('Sorting by energy done.')

    # writing the output to a .csv file, readable by excel
    with utilities.open_csv((prefix + '.csv').replace('..csv', '.csv'),
                            'w') as csvfile:  # open_csv chooses correct encoding based on whether python 2 or 3 is used
        write_csv = csv.writer(csvfile, delimiter=',', quotechar='|',
                               quoting=csv.QUOTE_MINIMAL)  # the str(u'') constructs are needed if unicode_literals are imported. See:https://github.com/jdunck/python-unicodecsv/issues/36
        write_csv.writerow(results_headers)
        for iterator in range(0, N_files):
            write_csv.writerow(results_summary.values[index_Esorted[iterator][0]])
        write_csv.writerow(['Structural information'])
        write_csv.writerow(
            ['A (Ang)', 'B (Ang)', 'C (Ang)', 'alpha', 'beta', 'gamma', 'cell_volume (Ang^3)', 'N_atoms', 'spg H-M'])
        write_csv.writerow(np.concatenate((
                                          resultsi_all[0].get_cell_lengths_and_angles(), [resultsi_all[0].get_volume()],
                                          [resultsi_all[0].get_number_of_atoms()], [spg])))
        write_csv.writerow(['Cartesian', 'Crystal', 'Axes'])
        for axes in resultsi_all[0].get_cell():
            write_csv.writerow(axes)
        write_csv.writerow(['Distances between final muon sites (Angstrom) '])
        for iterator in range(0, N_files):
            write_csv.writerow(muon_distances[iterator])  # the muon_distances array is already correctly sorted
        for iterator in range(0, N_files):
            write_csv.writerow([' '])
            for row_iterator in results_detailedpandas[index_Esorted[iterator][0]]:
                write_csv.writerow(row_iterator)
    if verbose in ['min', 'max']:
        print('... done.')
