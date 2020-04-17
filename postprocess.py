#!/usr/local/bin/python3.7
# postprocess.py - performs post-processing on DFT runs completed with pw.x
# (essentially just a wrapper for Franz's analysis.py)
# John Wilkinson 13/5/19

import sys
import os
sys.path.insert(0, os.path.abspath("/Users/johnny/Documents/University/DFT/qe_utils/"))
import Franz_DFT_scripts.analysis.analysis as analysis  # import Franz's analysis code


def main():
    # setup the parameters to use
    prefix = 'YF3.relax.mu'  # common prefix for all files
    postfix = '.pwo'  # common postfix (.out for now)
    spacegroup = 62
    makeCIFs = True
    muon_multiplicity_tol = 0.2 # maximum allowed distance to maintain multiplicity
    dist_cluster_tol = 0.5 # maximum distance between sites to be part of the same cluster

    analysis.QErel(prefix=prefix, postfix=postfix, set_spg=spacegroup, cif_flag=makeCIFs,
                   muon_multipl_prec=muon_multiplicity_tol, verbose='max', dist_cluster_prec=dist_cluster_tol,
                   supercell=[2,2,2])


if __name__=='__main__':
    main()
