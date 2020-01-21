# -*- coding: utf-8 -*-
"""
A few key utilities used in various other DFT+mu scripts.
Franz Lang
"""
from __future__ import print_function #python 2-3 compatible printing function
from ase import io #https://wiki.fysik.dtu.dk/ase/ase/atoms.html#list-of-all-methods
import numpy as np
import spglib, sys
import itertools, math

def progbar(i, i_max, bar_len=20):
    """A textual progress bar for the command line. Adapted from S. Sturniolo's Soprano scripts. Franz Lang 02/2018

    | Args:
    |   i (int): current progress index
    |   max_i (int): final progress index
    |   bar_len (Optional[int]): length in characters of the bar (no brackets)
    | Returns:
    |   bar (str): a progress bar formatted as requested

    """
    block = {True: 'X', False: ' '} #unicode \u2588 is the full block character
    perc = i/float(i_max)*bar_len
    bar = '[{0}]'.format(''.join([block[i < perc] for i in range(bar_len)]))

    return bar

def get_spg(file_or_structure, num_syms_flag=True,verbose=True):
    """A function to determine the spacegroup of a structure and the number of symmetry operations. Frang Lang 03/2018

    | Args:
    |   filename (str): name of file from which the structure will be loaded (e.g. .cif)
    |   num_syms_flag (bool): whether to calculate and return the number of symmetry operations
    |   verbose (bool): whether any progress should be printed
    | Returns:
    |   spacegroup (str): H-M symbol of the detected spacegroup
    |   num_syms (int): number of symmetry operations

    """
    #python 2-3 compatible check if a string (i.e. filename) was passed
    try:
        basestring
    except NameError:
        basestring = str
    if isinstance(file_or_structure, basestring): #if a string was passed to the function
        struct=io.read(file_or_structure) #read in the structure from the file
    else:#otherwise assume an ASE Atoms object was passed ot the program
        struct=file_or_structure #use the passed structure
    lattice = np.array(struct.get_cell().T, dtype='double', order='C')
    positions = np.array(struct.get_scaled_positions(), dtype='double', order='C')
    numbers = np.array(struct.get_atomic_numbers(), dtype='intc')
    cell = (lattice, positions, numbers)
    spg=spglib.get_spacegroup(cell)#spacegroup of system. H-M symbol followed by spacegroup number, e.g.: P2_1/c (14)
    if verbose:
        print('Spacegroup detected: ',spg)
    if num_syms_flag:
        num_syms=len(spglib.get_symmetry(cell)['rotations']) #number of symmetry operations in this spacegroup; equivalent to the maximum multiplicity.
        if verbose:
            print('Number of symmetry operations detected: ',str(num_syms))
        return spg,num_syms
    return spg

def file_strchange(filename, old_string, new_string,newfile_ext=None):
    """A function to change a string in a file to another string. Useful for changing an unconventional atom label
    such as 'Mu' to something ASE and spglib will understand, e.g. 'X' or 'XMu'. Frang Lang 03/2018

    | Args:
    |   filename (str): name of file to be changed
    |   old_string (str): string to be changed
    |   new_string (str): string that the 'old_string' will be replaced with
    |   newfile_ext (str): if set a copy of the original file will be made with the original file extension replaced with .'newfile_ext'.extension; if not set the original file will be changed.
    | Returns:
    |   nothing

    """
    # Safely read the input filename using 'with'
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            return
    if newfile_ext is not None:
        filename=filename.rsplit('.',1)
        filename=filename[0]+'.'+newfile_ext+filename[-1]
    # Safely write the changed content, if found in the file
    with open(filename, 'w') as f:
        s = s.replace(old_string, new_string)
        f.write(s)
    return

def open_csv(filename, mode='r'):
    """Open a csv file in proper mode depending on Python verion."""
    return(open(filename, mode=mode+'b') if sys.version_info[0] == 2 else
           open(filename, mode=mode, newline=''))
    
def wimda2gle(infile,plotlist,outfile=None,xcolumn=4,xmin=0,xmax=300,xtick=50,xstick=10,xtitle='T (K)',columns=1,expr=None):
    if outfile is None: #if user did not specify name of output file
        outfile=infile.rsplit('.',1)[0]+'.gle' #use input file name with extension changed to .gle
    Ngraphs=len(plotlist) #total number of graph panels to be plotted
    graphspercolumn=math.ceil(Ngraphs/columns) #number of graph panels per column (last column might have fewer)
    labelwidth=3.6 #width allowed for y-labels and y-title
    plotwidth=21/columns-labelwidth #width of the actual graphs (not including labels or titles as fullsize setting is used)
    plotheight=(29.7-3.5)/graphspercolumn #height of the actual graphs
    #preamble setting various paramters for the gle plot
    glepreamble=('size 21 29.7\n\
ww=%s\n'%(plotwidth)+\
'hi=%s\n'%(plotheight)+\
'xsep=%s\n'%(labelwidth)+\
'set font texcmr\n\
set texlabels 1\n\
set texscale scale\n\
ht=1.0\n\
htkey=1.0\n\
set hei ht\n\
ms=0.3\n\
lw=0.05\n\
xmin=%s\n'%(xmin)+\
'xmax=%s\n'%(xmax)+\
'xtik=%s\n'%(xtick)+\
'xstik=%s\n'%(xstick)+\
'color0=\'blue\'\n\
color1=\'red\'\n\
color2=\'green\'\n\
color3=\'orange\'\n\
color4=\'black\'\n\
begin translate 2 2.1\n\
')
    gleend='end translate' #final line, needed for ending the translate environment
    graphs='' #initialise string that will contain all the graphs
    for i in range(0,Ngraphs): #iterate through all the graphs
        columns=list(itertools.chain.from_iterable(plotlist[i][1])) #list of the data file columns to plotted in the current graph
        xmove=math.floor(i/graphspercolumn) #current column/x-position of current panel
        ymove=graphspercolumn-i%graphspercolumn-1 #y-position of current panel via modulo operator
        graphs=graphs+\
'amove 1+%s*(ww+xsep) %s*hi\n'%(xmove,ymove)+\
'begin graph \n\
size ww hi\n\
fullsize\n\
xaxis min xmin max xmax dticks xtik dsubticks xstik\n\
ytitle "%s" dist 0.1\n'%(plotlist[i][2])
        if ymove!=0 and i!=Ngraphs-1: #x-labels off for graphs that are not at the bottom of their column
            graphs=graphs+'xlabels off\n'
        if ymove!=graphspercolumn-1: #last y-label off for graphs that are not at the top of their column
            graphs=graphs+'yaxis nolast\n'
        if ymove==0 or i==Ngraphs-1: #x-title for graphs that are at the bottom of their column
            graphs=graphs+'xtitle "%s" dist 0.1 hei ht\n'%(xtitle)
        graphs=graphs+'data %s '%(plotlist[i][0]) #data file to be read in
        for j,column in enumerate(columns):
            graphs=graphs+'d%s=c%s,c%s '%(j+1,xcolumn,column) #data columns to be read in; x-column=4 is default for T values in wimda fit tables 
        graphs=graphs+'\n'
        counter=1 #counter variable making sure that correct quantities are plotted if no error data is given
        for j,colpairs in enumerate(plotlist[i][1]):
            if len(colpairs)==2:
                graphs=graphs+'d%s err d%s marker fcircle color %s lwidth lw msize ms\n'%(counter,counter+1,'color'+str(j)) #plot data including errors
                counter=counter+2
            elif len(colpairs)==1:
                graphs=graphs+'d%s marker fcircle color %s lwidth lw msize ms\n'%(j,'color'+str(j)) #plot data without associated errors
                counter=counter+1
        graphs=graphs+'begin key\n\
nobox\ndist 0.2\ncoldist 0.2\npos br\nmargins 0 0\njustify br\noffset ww-0.2 0.3\n' #write the preamble for the graph key
        for j,keys in enumerate(plotlist[i][3]):
            if j==0:
                graphs=graphs+'marker fcircle msize ms color color0 text %s\n'%('"'+keys+'"') #write the first key
            else:
                graphs=graphs+'separator\nmarker fcircle msize ms color color%s text %s\n'%(j,'"'+keys+'"') #write the other key entries
        graphs=graphs+'end key\nend graph\n' #end key and graph
    if expr is not None:
        graphs=graphs+'amove 1 %s*hi+0.4\n'%(graphspercolumn)+'tex %s\n'%(expr)
    try: #open the file and write the gle-file
        wfile=open(outfile,'w')
        wfile.write(glepreamble)
        wfile.write(graphs)
        wfile.write(gleend)
        wfile.close()
    except IOError:
        print('Could not open file %s'%(outfile))