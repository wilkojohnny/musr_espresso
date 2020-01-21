# Soprano - a library to crack crystals! by Simone Sturniolo
# Copyright (C) 2016 - Science and Technology Facility Council

# Soprano is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Soprano is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
selection.py

Contains the definition of an AtomSelection class,
namely a group of selected atoms for a given structure,
and methods to build it.
"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import hashlib
import warnings
import numpy as np

from soprano.utils import minimum_supcell, supcell_gridgen


# This decorator applies to all operators providing some basic checks
def _operator_checks(opfunc):

    def decorated_opfunc(self, other):
        if not isinstance(other, AtomSelection):
            raise TypeError('AtomSelection does not support operations with'
                            ' different types')

        if self._auth is not None and other._auth is not None:
            # Check compatibility
            if self._auth != other._auth:
                raise ValueError('Selections come from different systems')

        return opfunc(self, other)

    return decorated_opfunc


class AtomSelection(object):

    """AtomSelection object.

    An AtomSelection represents a group of atoms from an ASE Atoms object.
    It keeps track of them and can be used to perform operations on them
    (for example geometrical transformation or extraction of specific
    properties).
    It does not keep track of the original Atoms object it's been created
    from, but can be "authenticated" to verify that it is indeed operating
    consistently on the same structure. It also provides a series of static
    methods to build selections with various criteria.

    """

    def __init__(self, atoms, sel_indices, authenticate=True):
        """Initialize the AtomSelection.

        | Args:
        |   atoms (ase.Atoms): the atoms object on which the selection is
        |                      applied
        |   sel_indices (list[int]): the list of indices of the atoms that
        |                            are to be selected
        |   authenticate (Optional[bool]): whether to use hashing to confirm
        |                                  the identity of the atoms object
        |                                  we're operating with

        """

        # A quick check: are the indices actually contained in the Atoms?
        if len(sel_indices) > 0:
            if (min(sel_indices) < 0 or
                    max(sel_indices) >= atoms.get_number_of_atoms()):
                raise ValueError('Invalid indices for given Atoms object')

        self._indices = np.array(sel_indices)

        if authenticate:
            # Create an hash for certification
            self._auth = self._hash(atoms)
        else:
            self._auth = None

        self._arrays = {}

    @property
    def indices(self):
        return self._indices

    def _hash(self, atoms):
        """A function to create an identifying hash for a given Atoms system.
        This is used later to check that the system is indeed unchanged when
        the Selection is reused.

        While changes in positions or cell don't invalidate the Selection,
        changes in composition potentially do (indices of atoms can change).
        """

        h = hashlib.md5()
        h.update(''.join(atoms.get_chemical_symbols()).encode())

        return h.hexdigest()

    def set_array(self, name, array):
        """Save an array of given name containing arbitraty information
        tied to the selected atoms.
        This must match the length of the selection and will be passed on to
        any Atoms objects created with .subset.

        | Args:
        |   name (str): name of the array to be set or created
        |   array (np.ndarray): array of data to be saved

        """

        # First a check
        if len(array) != len(self):
            raise ValueError("Invalid array passed to set_array")

        self._arrays[name] = np.array(array)

    def get_array(self, name):
        """Retrieve a previously stored data array.

        | Args:
        |   name (str): name of the array to be set or created

        | Returns:
        |   array (np.ndarray): array of data to be saved

        """

        # If the name isn't right just let the KeyError happen
        return self._arrays[name]

    def validate(self, atoms):
        """Check that the given Atoms object validates with this selection."""
        if self._auth is None:
            warnings.warn('WARNING'
                          ' - this selection does not support validation')
            return True
        else:
            return self._hash(atoms) == self._auth

    def subset(self, atoms):
        """Generate an Atoms object containing only the selected atoms."""

        if not self.validate(atoms):
            raise ValueError(
                'Given Atoms object does not match this selection')

        subset = atoms.copy()
        not_sel = list(set(range(atoms.get_number_of_atoms())) -
                       set(self._indices))
        del subset[not_sel]

        # Now the arrays
        for k in self._arrays:
            # The array needs to be sorted with the indices
            # in order to match the order of the atoms in the Atoms object
            _, arr = zip(*sorted(zip(self._indices, self._arrays[k]),
                                 key=lambda t: t[0]))
            subset.new_array(k, np.array(arr))

        return subset

    def __getitem__(self, indices):
        """Slicing: take only part of a selection"""

        if type(indices) is int:
            # Special case, a single element!
            indices = slice(indices, indices+1)

        try:
            newsel = self._indices[indices]
        except TypeError:
            newsel = [self._indices[i] for i in indices]

        sliced = copy.deepcopy(self)
        sliced._indices = newsel
        sliced._arrays = {k: a[indices] for k, a in self._arrays.iteritems()}

        return sliced

    def __iter__(self):
        return [self[i:i+1] for i in range(len(self))].__iter__()

    # Overloading operators to allow sum, subtraction and product of selections
    @_operator_checks
    def __add__(self, other):
        """Sum: join selections"""

        # Join
        ans = copy.deepcopy(self)
        ans._indices = np.array(list(set(self.indices).union(other.indices)))
        # For the arrays:
        # only join the ones present in BOTH selections
        common_k = set(self._arrays.keys()
                       ).intersection(set(other._arrays.keys()))
        ans._arrays = {}
        for k in common_k:
            ans._arrays[k] = np.concatenate((self._arrays[k],
                                             other._arrays[k]))

        return ans

    @_operator_checks
    def __sub__(self, other):

        # Difference
        ans = copy.deepcopy(self)
        ans._indices = np.array(list(set(self.indices)-set(other.indices)))
        # For the arrays:
        # keep them but remove the removed indices
        arr_i = [np.where(self.indices == i)[0][0] for i in ans._indices]
        for k in ans._arrays:
            ans._arrays[k] = ans._arrays[k][arr_i]

        return ans

    @_operator_checks
    def __mul__(self, other):

        # Intersection
        ans = copy.deepcopy(self)
        ans._indices = np.array(list(set(self.indices)
                                     .intersection(other.indices)))
        # For the arrays:
        # keep the ones present in either selection,
        # but only the relevant indices of course,
        # and remove if conflicting!
        all_k = set(self._arrays.keys()).union(set(other._arrays.keys()))
        arr1_i = [np.where(self.indices == i)[0][0] for i in ans._indices]
        arr2_i = [np.where(other.indices == i)[0][0] for i in ans._indices]
        ans._arrays = {}
        for k in all_k:
            try:
                arr1 = self._arrays[k][arr1_i]
            except KeyError:
                arr1 = None
            try:
                arr2 = other._arrays[k][arr2_i]
            except KeyError:
                arr2 = None

            if arr1 is not None and arr2 is not None:
                # Do they conflict?
                if not np.all(arr1 == arr2):
                    print(('WARNING - conflicting arrays of name {0} found'
                           ' will be removed during intersection'
                           ' operation').format(k))
                    continue

            ans._arrays[k] = arr1 if arr1 is not None else arr2

        return ans

    def __len__(self):
        return len(self._indices)

    def __contains__(self, item):
        return item in self._indices

    @staticmethod
    def all(atoms):
        """Generate a selection for the given Atoms object of all atoms.

        | Args:
        |   atoms (ase.Atoms): Atoms object on which to perform selection

        | Returns:
        |   selection (AtomSelection)

        """

        return AtomSelection(atoms, range(len(atoms)))

    @staticmethod
    def from_element(atoms, element):
        """Generate a selection for the given Atoms object of all atoms of a
        specific element.

        | Args:
        |   atoms (ase.Atoms): Atoms object on which to perform selection
        |   element (str): symbol of the element to select

        | Returns:
        |   selection (AtomSelection)

        """

        sel_i = np.where(np.array(atoms.get_chemical_symbols()) == element)[0]

        return AtomSelection(atoms, sel_i)

    @staticmethod
    def from_box(atoms, abc0, abc1, periodic=False, scaled=False):
        """Generate a selection for the given Atoms object of all atoms within
        a given box volume.

        | Args:
        |   atoms (ase.Atoms): Atoms object on which to perform selection
        |   abc0 ([float, float, float]): bottom corner of box
        |   abc1 ([float, float, float]): top corner of box
        |   periodic (Optional[bool]): if True, include periodic copies of the
        |                              atoms
        |   scaled (Optional[bool]): if True, consider scaled (fractional)
        |                            coordinates instead of absolute ones

        | Returns:
        |   selection (AtomSelection)

        """

        if scaled:
            pos = atoms.get_scaled_positions()
        else:
            pos = atoms.get_positions()
        # Do we need periodic copies?
        if periodic and any(atoms.get_pbc()):
            # Get the range
            max_r = np.linalg.norm(np.array(abc1)-abc0)
            scell_shape = minimum_supcell(max_r, latt_cart=atoms.get_cell(),
                                          pbc=atoms.get_pbc())
            grid_frac, grid = supcell_gridgen(atoms.get_cell(), scell_shape)
            if scaled:
                pos = (pos[:, None, :]+grid_frac[None, :, :])
            else:
                pos = (pos[:, None, :]+grid[None, :, :])

        where_i = np.where(np.all(pos > abc0, axis=-1) &
                           np.all(pos < abc1, axis=-1))[:2]

        sel_i = where_i[0]

        sel = AtomSelection(atoms, sel_i)
        if periodic:
            sel.set_array('cell_indices', grid_frac[where_i[1]])

        return sel

    @staticmethod
    def from_sphere(atoms, center, r, periodic=False, scaled=False):
        """Generate a selection for the given Atoms object of all atoms within
        a given spherical volume.

        | Args:
        |   atoms (ase.Atoms): Atoms object on which to perform selection
        |   center ([float, float, float]): center of the sphere
        |   r (float): radius of the sphere
        |   periodic (Optional[bool]): if True, include periodic copies of the
        |                              atoms
        |   scaled (Optional[bool]): if True, consider scaled (fractional)
        |                            coordinates instead of absolute ones

        | Returns:
        |   selection (AtomSelection)

        """

        if scaled:
            pos = atoms.get_scaled_positions()
        else:
            pos = atoms.get_positions()
        # Do we need periodic copies?
        if periodic and any(atoms.get_pbc()):
            # Get the range
            r_bounds = minimum_supcell(r, latt_cart=atoms.get_cell(),
                                       pbc=atoms.get_pbc())
            grid_frac, grid = supcell_gridgen(atoms.get_cell(), r_bounds)
            if scaled:
                pos = (pos[:, None, :]+grid_frac[None, :, :])
            else:
                pos = (pos[:, None, :]+grid[None, :, :])

        where_i = np.where(np.linalg.norm(pos-center, axis=-1) <= r)

        sel_i = where_i[0]

        sel = AtomSelection(atoms, sel_i)
        if periodic:
            sel.set_array('cell_indices', grid_frac[where_i[1]])

        return sel
