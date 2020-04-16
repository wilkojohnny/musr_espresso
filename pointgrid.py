# -*- coding: utf-8 -*-
from __future__ import print_function #python 2-3 compatible printing function
import timeit, sys
import numpy as np
from ase import io, Atoms, Atom, spacegroup, visualize  # https://wiki.fysik.dtu.dk/ase/ase/atoms.html#list-of-all-methods
import utilities
import copy


def point_grid(input_structure,xgrid=1,ygrid=1,zgrid=1,random=None,rand_factor=1.1,round_to=None,set_spg=None,set_spgsetting=None,set_numsyms=None,prec_atom=0.5,prec_grid=0.1,verbose='min',export_flag=False,export_file=None,number_format='%.3f',out_grid=True,out_struct=False,export_cif=False):
    """A function to create a grid points in a given spacegroup/structure that are both a minimum distance away from all atom
    in the structure and from each other. Frang Lang 03/2018

    | Args:
    |   filename (str): name of file from which the structure will be loaded (e.g. .cif) (JW-REMOVED)
    |   input_structure: either ase structure or a filename of a structure file. (JW)
    |   xgrid, ygrid, zgrid (int): number of grid points in each cartesian direction
    |   random (int): if set a starting grid of a randomly generated set of points will be used, where `random` is the number of starting points
    |   rand_factor (float>1): factor by which the number of randomly chosen points is multiplied. This can be used to create exactly 'random' number of points: the program creates int(random*rand_factor) of points and then deletes those due to symmetries, too close to atoms etc and if the remaining number of points > 'random' it truncates the list to exactly 'random' number of points.
    |   round_to (int): number of decimals to which the generated points are rounded (might be needed to avoid making points the same when exporting them in a format that will truncate the decimals)
    |   set_spg (str or int): set spacegroup explicitly by either H-M symbol or number in Int. Tables; if not set the function will determine it from the structure file
    |   set_spgsetting (int): use a non-conventional setting for a spacegroup; if not set default value of 1 will be used
    |   set_numsyms (int): the number of symmetry operations; if not set the function will determine it from the structure file
    |   prec_atom (float): minimum allowed distance (in Angstrom) between grid points and their nearest atoms
    |   prec_grid (float): minimum allowed distance (in Angstrom) between grid points
    |   verbose ('min','max'): whether any progress should be printed
    |   export_flag (bool): flag whether the the resultant grid should be exported or not
    |   export_file (str): name of file to export the results to; if not set will use 'filename' with the extension replaced by 'xgridxygridxzgridgrid.txt'
    |   number_format (str): formatting of the numbers when they are exported to a text file
    |   export_cif (bool): flag whether to export the grid of points as a cif file (name is `filename` with extension replaced by `xgridxygridxzgridgrid.cif`)
    | Returns:
    |   if out_grid=True: grid_survivors: the low symmetry grid points that satisfy a minimum distance to all atoms and all other grid points
    |   if out_struct=True: structure read in from `filename` as an ASE atoms object

    """
    # JW - change filename argument to input_structure, so that premade structures can be given in the argument
    if isinstance(input_structure, str):
        struct = io.read(input_structure)  # read in structure from file
        filename = input_structure
    else:
        struct = input_structure
        filename = 'structure'

    N_max_approx=(struct.get_volume()-4*3.1415/3*prec_atom**3*struct.get_number_of_atoms())/(4*3.1415/3*prec_grid**3) #approximate number of maximum grid points
    if verbose in ['min', 'max']:
        print('Approximate maximum number of grid points based on prec_atom and prec_grid: ',str(int(N_max_approx)))

    # determine spacegroup and setting (dependent on whether user has set any by hand)
    spg, num_syms=utilities.get_spg(struct, num_syms_flag=True, verbose=False)
    spg = spg.split()[0] #H-M of spacegroup only. Element [1] would be number in international tables
    if set_spg is None:
        if verbose in ['min','max']:
            print('Spacegroup detected from file: ',str(spg))
        if set_numsyms is None:
            if verbose in ['min','max']:
                print('Number of symmetry operations from file: ',str(num_syms))
        else:
            num_syms=set_numsyms
            if verbose in ['min','max']:
                print('Number of symmetry operations set: ',str(num_syms))
    else:
        spg=set_spg
        if verbose in ['min','max']:
            print('Spacegroup set by user: ',str(spg))
        if set_numsyms is None:
            if verbose in ['min','max']:
                print('Number of symmetry operations from file: ',str(num_syms))
        else:
            num_syms=set_numsyms
            if verbose in ['min','max']:
                print('Number of symmetry operations set: ',str(num_syms))
    if set_spgsetting is None:
        spgsetting=1
        if verbose in ['min','max']:
            print('Default setting: 1')
    else:
        spgsetting=set_spgsetting
        if verbose in ['min','max']:
            print('User setting: ',str(spgsetting))

    # create the full grid of points depending on whether a regular or a random grid is wanted
    # regular grid: create xgrid*ygrid*zgrid regular grid of points
    if random is None:
        N_points = xgrid*ygrid*zgrid #total number of starting gridpoints
        grid_full = np.zeros((N_points,3),dtype=float) #initialise array to hold all gridpoints
        counter=0
        for xi in range(0,xgrid):
            for yi in range(0,ygrid):
                for zi in range(0,zgrid):
                    grid_full[counter,:] = [1.*xi / (xgrid), 1.*yi / (ygrid), 1.*zi/(zgrid)] #create uniform grid of points
                    counter = counter+1
    else: # random grid: grid of random points
        N_points=int(random*rand_factor) # total number of starting gridpoints. Use 'rand_factor' to overestimate, so that there is a real chance of at least 'random' number of points being left at the end.
        grid_full = np.random.rand(N_points, 3) # create random grid of points
    # round the points of the starting grid if set so by user
    if round_to is not None:
        grid_full=np.round(grid_full,round_to)
    # create a list which functions as a selector for which points will be used
    grid_selector=np.array([False for i in range(0,N_points)],dtype=bool)
    #start a timer
    start_time = timeit.default_timer()

    # read in the structure again to create an independent copy
    if isinstance(input_structure, str):
        struct_dump = io.read(input_structure)
    else:
        struct_dump = copy.deepcopy(input_structure)

    # look at each gridpoint in turn
    point_crystal_all=None

    #check multiplicity of each grid point
    for counter,grid_point in enumerate(grid_full):
        #write a progress bar every 10 steps unless `max` verbose output is wanted (as it will interfere with the other printed output)
        if verbose in ['min']:
            if counter==0 or counter%10==9 or counter==N_points-1:
                sys.stdout.write("\rPoint: "+str(counter+1).rjust(len(str(N_points)))+'/'+str(N_points)+' {0}'.format(utilities.progbar(counter+1,N_points))+' t='+str(np.round(timeit.default_timer() - start_time,2))+'secs')

        point_crystal=spacegroup.crystal(symbols='X',basis=grid_point,spacegroup=spg,setting=spgsetting,cell=struct.get_cell()) #create crystal with grid point as basis
        point_multiplicity=point_crystal.get_number_of_atoms() #multiplicity of the grid point
        if verbose in ['max']:
            print('Multiplicity of point ',str(grid_point),': ',point_multiplicity)
        if True: #point_multiplicity==num_syms: #JW - remove the multiplicity thing for now... #if grid point has lowest symmetry possible we need to check if it's too close to any of the atoms
            if verbose in ['max']:
                print('low sym. point: ',str(grid_point))

            #check distance of grid point to atoms of system
            point_atoms=Atoms(symbols=['X'],scaled_positions=[grid_point],cell=struct.get_cell()) #create an Atoms object for the grid point
            struct_dump.extend(point_atoms) #add the atoms object created from the grid point
            point_atomdistances=struct_dump.get_distances(-1,range(0,struct_dump.get_number_of_atoms()-1),mic=True) #distances of all atoms and already kept grid points from the grid point at hand
            struct_dump.pop() #remove the grid point from the original structure
            #print 'smallest atom dist.: '+str(np.min(point_atomdistances))
            if np.min(point_atomdistances)>prec_atom: #if the grid point is further than prec_atom (in Angstrom) from all atoms check further: its distance to already selected grid points
                if point_crystal_all is None: #if it's the first low symmetry grid point initialise a lattice of grid points and set corresponding selector to True
                    point_crystal_all=point_crystal
                    grid_selector[counter]=True
                    if verbose in ['max']:
                        print('1st grid point: ',str(grid_point))

                #check distance of grid point to other grid points
                else:
                    point_crystal_all.extend(point_crystal) #add the crystal created from the current grid point to the existing lattice of already kept grid points
                    point_pointdistances=point_crystal_all.get_distances(-num_syms,range(0,point_crystal_all.get_number_of_atoms()-num_syms),mic=True) #get distances to already kept grid points
                    if verbose in ['max']:
                        print('smallest grid point dist.: ',str(np.min(point_pointdistances)))
                    if np.min(point_pointdistances)>prec_grid: #if the grid point is further than prec_grid (in Angstrom) from all already kept grid points then keep it as well
                        grid_selector[counter]=True
                        if verbose in ['max']:
                            print('another grid point: ',str(grid_point))
                    else:
                        for i in range(0,num_syms):
                            point_crystal_all.pop() #if current grid point is closer than prec_grid to already kept grid points then delete it from grid point lattice
    grid_survivors=grid_full[grid_selector] #select the 'surviving' grid points, which are: i)on a lowest symmetry point, ii) at least a distance prec_atom (prec_grid) away from all other atoms (grid points)
    if random is not None and random<len(grid_survivors):#double-check if the number of 'surviving' points is larger than the requested number of random points. May happen if 'rand_factor'>1.
        grid_survivors=grid_survivors[:random]
    if verbose in ['min','max']:
        print('\n',str(xgrid),'x',str(ygrid),'x',str(zgrid),' grid for spacegroup ',str(spg),' (setting ',str(spgsetting),') gave ',str(len(grid_survivors)),' grid points.')
    if export_flag: #if text file with grid points should be exported
        if export_file is None: #if no explicit name for the export file is set by user
            if random is None: #if a regular grid was created
                export_filename=filename.rsplit('.',1)[0]+'_'+str(xgrid)+'x'+str(ygrid)+'x'+str(zgrid)+'grid.txt'
            else: #if a grid of random points was created
                export_filename=filename.rsplit('.',1)[0]+'_'+str(len(grid_survivors))+'randgrid.txt'
        else:
            export_filename=export_file
        if verbose in ['min','max']:
            print('Exporting grid points to ',export_filename,'.')
        print(grid_survivors)
        np.savetxt(export_filename,grid_survivors,fmt=number_format)
    if export_cif: #if grid points are supposed to be exported as a cif-file
            grid_cifexport=Atoms(symbols=['X' for i in range(0,len(grid_survivors))],scaled_positions=grid_survivors,cell=struct.get_cell())
            if random is None: #if a regular grid was created
                export_cifname=filename.rsplit('.',1)[0]+'_'+str(xgrid)+'x'+str(ygrid)+'x'+str(zgrid)+'grid.cif'
            else: #if a grid of random points was created
                export_cifname=filename.rsplit('.',1)[0]+'_'+str(len(grid_survivors))+'randgrid.cif'
            io.write(export_cifname,grid_cifexport)
            print('Exporting grid points to ',export_cifname,'.')
    if out_grid:
        if out_struct:
            return grid_survivors,struct
        else:
            return grid_survivors
    else:
        if out_struct:
            return struct
        else:
            return


def get_internuclei_positions(crystal_cell: Atoms, nn_nuclei: list = None, max_distance: float = 3.0, spg=None, prec_grid: float = 0.1) -> list:
    """
    Get a set of muon positions which are in between the nuclei with labels nn_nuclei. John Wilkinson, 16/4/2020
    :param crystal_cell: cell of crystal to put the nuclei into
    :param nn_nuclei: list of strings of the labels of the nn nuclei. Default is ['F'].
    :param max_distance: maximum possible nn nuclei distance from each other
    :param spg: space group of crystal
    :param prec_grid: maximum distance between two muons during the symmetry analysis to be treated as the same
    :return: list of SCALED positions of muons which are in between the nn_nuclei nuclei
    """

    if nn_nuclei is None:
        nn_nuclei = ['F']

    potential_positions = Atoms(cell=crystal_cell.get_cell())

    # calculate all the cell distances initially (for speed)
    cell_distances = crystal_cell.get_all_distances(mic=True)

    # for each atom in crystal_cell
    for i_atom in range(0, len(crystal_cell)):
        # is this atom one of the ones we are interested in?
        if crystal_cell[i_atom].symbol in nn_nuclei:
            # get all the distances between this nuclei and the others
            i_atom_distances = cell_distances[i_atom]
            for j_atom in range(i_atom, len(i_atom_distances)):
                # for each atom, if it has the correct symbol, *and* the distance is acceptable:
                if crystal_cell[j_atom].symbol in nn_nuclei and max_distance > i_atom_distances[j_atom] > 0:
                    this_potential_position = (crystal_cell[i_atom].position + crystal_cell[j_atom].position) / 2
                    potential_positions.append(atom=Atom('H', this_potential_position, tag=i_atom+j_atom*len(crystal_cell)))

    # now seek and destroy duplicates
    # for each muon position
    potential_positions_scaled = potential_positions.get_scaled_positions()
    remove_list = []
    for i_muon_position, muon_position in enumerate(potential_positions_scaled):
        if i_muon_position not in remove_list:
            sym_muon_positions_i = spacegroup.crystal('H', basis=muon_position, cell=potential_positions.get_cell(),
                                                      spacegroup=spg)

            for j_muon_position in range(i_muon_position+1, len(potential_positions)):
                if j_muon_position not in remove_list:  # i.e if this one isn't already marked for removal
                    sym_muon_positions_j = spacegroup.crystal('H', basis=potential_positions_scaled[j_muon_position],
                                                              cell=potential_positions.get_cell(), spacegroup=spg)

                    # now try to merge the two -- and if any of the distances are less than prec_grid, then rm j later
                    sym_muon_positions_ij = copy.deepcopy(sym_muon_positions_i)
                    for atom in sym_muon_positions_j:
                        sym_muon_positions_ij.append(atom)

                    dist = sym_muon_positions_ij.get_all_distances(mic=True)
                    # this next line, in the ugliest way possible, finds the minimum distance above 0 (0 is
                    min_dist = [min(dist[i][j] for i in range(0, len(dist[:])) for j in range(0, len(dist[:])) if i > j)]
                    if min_dist[0] < prec_grid:
                        print(i_muon_position)
                        print(j_muon_position)
                        remove_list.append(j_muon_position)

    # now destroy duplicates
    remove_list.sort(reverse=True)
    print(remove_list)
    for remove_item in remove_list:
        potential_positions.pop(remove_item)
        potential_positions_scaled = np.delete(potential_positions_scaled, remove_item, axis=0)

    # visualize.view(potential_positions)

    return potential_positions_scaled