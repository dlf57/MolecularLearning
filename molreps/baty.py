'''
Function for reading common molecule files and 
creating a Typed Bond Angle representations.
'''

# Author: Dakota Folmsbee <dfolmsbee@gmail.com>
# License: GPLv2

import pybel
import glob
import copy
import numpy as np
from itertools import chain
from math import sqrt
from math import radians

# don't need to reimport openbabel
ob = pybel.ob


def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()


encode = {'Br': '350', 'C+': '064', 'C1': '061', 'C2': '062',
          'C3': '063', 'Cac': '062', 'Car': '062',
          'Cl': '170', 'F': '090', 'H': '010', 'HO': '010',
          'I': '530', 'N1': '071', 'N2': '072',
          'N3': '073', 'N3+': '073', 'Nam': '073',
          'Nar': '072', 'Ng+': '072', 'Nox': '072',
          'Npl': '072', 'Ntr': '071', 'O-': '085',
          'O.co2': '082', 'O2': '082', 'O3': '083',
          'P': '150', 'S2': '162', 'S3': '163',
          'Sac': '162', 'Se': '340', 'So2': '162',
          'Sox': '162'}
'''Atom encoding
Each atom encoding is represented by a string of numbers
consisting of an element encoding and a hybrid encoding.
Element encoding is done based off of the periodic table
where the element that will be encoded is encoded to the
corresponding element number on the periodic table.
Type encoding is done based off of the type of atom it is.
 _________________________
| Element (XX) | Type (X) |
Ex.
 ____________       ________
| Carbon | 3 | --> | 06 | 3 | --> 063
'''


def baty(mol_file, bags, bag_sizes):
    '''
    Function to iterate over the molecules and fill the
    bags with the appropriate values in order to construct
    the representations for machine learning.
    Paramters
    ---------
    mol_file: file
        molecule file for reading in coordinates
    bags: dict
        dict of all bags for the dataset
    bag_sizes: dict
        dict of size of the largest bags in the dataset
    Returns
    -------
    baty: vector
        vector of all bonds in the molecule
    '''
    # copy bags dict to ensure it does not get edited
    bag_set = copy.deepcopy(bags)

    # read molfile into openbabel
    filetype = mol_file.split('.')[1]
    mol = next(pybel.readfile(filetype, mol_file))

    # iterate through all bonds
    for bond in ob.OBMolBondIter(mol.OBMol):
        begin = atomType(mol, bond.GetBeginAtomIdx())
        end = atomType(mol, bond.GetEndAtomIdx())

        if (end < begin):
            # swap them for lexographic order
            begin, end = end, begin

        bond_type = encode[begin] + encode[end]
        bond_length = bond.GetLength()
        bag_set[bond_type].append(np.float16(bond_length))

    # iterate through all angles
    for angle in ob.OBMolAngleIter(mol.OBMol):
        a = (angle[0] + 1)
        b = mol.OBMol.GetAtom(angle[1] + 1)
        c = (angle[2] + 1)

        aType = atomType(mol, a)
        cType = atomType(mol, c)
        if (cType < aType):
            # swap them for lexographic order
            aType, cType = cType, aType

        angle_type = encode[aType] + encode[b.GetType()] + encode[cType]
        angle_angle = radians(abs(b.GetAngle(a, c)))
        bag_set[angle_type].append(np.float16(angle_angle))

    # sort bags by magnitude, pad, concactenate
    baty = []
    bag_keys = list(bag_set.keys())
    for i in range(len(bag_keys)):
        size = bag_sizes[bag_keys[i]] + 1
        baglen = len(bag_set[bag_keys[i]])
        if baglen > (size - 1):
            raise Exception(
                '{}-bag size is too small. '
                'Increase size to {}.'.format(bag_keys[i], baglen))
        pad = size - baglen
        bag_set[bag_keys[i]] = sorted(bag_set[bag_keys[i]], reverse=True)
        bag_set[bag_keys[i]].extend([0.] * pad)
        baty.append(bag_set[bag_keys[i]])

    # flatten bty into one list and store as a np.array
    baty = np.array(list(chain.from_iterable(baty)))

    return baty