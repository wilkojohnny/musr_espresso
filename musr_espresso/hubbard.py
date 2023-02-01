"""
Hubbard.py -- deals with hubbard parameters in pw.x

Created by John Wilkinson, 1/2/2023
"""

import math

class Hubbard(object):

    def __init__(self, atoms, hubbard_type='ortho-atomic'):
        # initially set the hubbard parameters to something low
        self.hubbard_type = hubbard_type
        self.U_atoms = {}
        atomic_numbers = atoms.get_atomic_numbers()
        for i_Z, Z in enumerate(atomic_numbers):
            # if Z>20, then will have at least d-electrons so likely to need a U
            if Z > 20 and atoms[i_Z] not in self.U_atoms:
                manifold_guess = math.floor((Z-21)/18) + 3
                self.U_atoms.update({atoms[i_Z]: {'{:d}d'.format(manifold_guess) : 1e-4}})

    def qe_7_hubbard(self, pp_names):
        """
        returns the hubbard parameters in a format for qe7
        """

        hubbard_string = 'HUBBARD {}\n'.format(self.hubbard_type)
        for atom in pp_names:
            # strip the number at the end of the pp_name to get the atom name
            atom_realname = ''.join(i for i in atom if not i.isdigit())
            try:
                pars = self.U_atoms[atom_realname]
                for manifold, Uval in pars.items():
                    hubbard_string += 'U {}-{} {} \n'.format(atom, manifold, Uval)
            except KeyError:
                print('No Hubbard U given for ' + atom_realname)

        return hubbard_string
