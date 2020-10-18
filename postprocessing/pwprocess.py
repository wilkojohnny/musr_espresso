#!/usr/local/bin/python3.7
import argparse
import os
import sys

from musr_espresso import analysis  # import Franz's analysis code


def main():
    # get the arguments
    parser = argparse.ArgumentParser(description='pw.x output file processor to find the muon site')
    parser.add_argument('-sg', '--spacegroup', nargs=1, default='1')
    parser.add_argument('-sc', '--supercell', nargs=3, type=int, required=False, default=[2, 2, 2],
                        help='Size of supercell used in the calculation, in terms of the number of unit cells'
                             ' (e.g 2 2 2)')
    parser.add_argument('-mtol', '--multiplicitytol', nargs=1, default=0.1, required=False, type=float,
                        help='Multiplicity tolerance -- maximum allowed distance (scaled coordinates) to maintain '
                             'multiplicity')
    parser.add_argument('-dctol', '--dclustertol', nargs=1, default=0.1, required=False, type=float,
                        help='Distance cluster tolerance -- maximum allowed distance (Angstroms) for muons to be in the'
                             'same cluster')
    parser.add_argument('pwofiles', nargs='*', type=str, help='pw.x output files to analyse')

    arguments = parser.parse_args()

    spacegroup = arguments.spacegroup
    if isinstance(spacegroup, list):
        spacegroup = spacegroup[0]
    if spacegroup.isdigit():
        spacegroup = int(spacegroup)
    supercell = arguments.supercell
    mtol = arguments.multiplicitytol
    dtol = arguments.dclustertol
    pwo_files = arguments.pwofiles

    try:
        assert len(pwo_files) > 0
    except AssertionError:
        print('No pw.x output files given. Bye. 😫')
        return 1

    analysis.QErel(pwo_files=pwo_files, set_spg=spacegroup, cif_flag=True,
                   muon_multipl_prec=mtol, verbose='max', dist_cluster_prec=dtol,
                   supercell=supercell)


if __name__ == '__main__':
    main()
