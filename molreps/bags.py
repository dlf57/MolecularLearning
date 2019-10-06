'''
Bag making class for bag of features style representations

TODO:
    - Add nonbonding of each
'''

# Author: Dakota Folmsbee <dfolmsbee@gmail.com>
# License: GPLv2

import copy
import glob
import pybel
from collections import OrderedDict

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


class BagMaker:
    """
    Class to make bags for bag representations
    Attributes
    ----------
    rep_str : str
        name of representation (ie. 'BATTY')
    dataset : path
        path to all molecules in the dataset
    """
    __accepted_reps = ['BTY', 'BTYNB', 'BATY', 'BATYNB', 'BATTY', 'BATTYNB']

    def __init__(self, rep_str=None, dataset=None):
        if (rep_str and dataset) is not None:
            self.rep(rep_str, dataset)
        return None

    def rep(self, rep_str, dataset):
        if rep_str == 'BTY':
            self.bty(dataset)
        elif rep_str == 'BATY':
            self.baty(dataset)
        elif rep_str == 'BATTY':
            self.batty(dataset)
        else:
            accept_reps = str(BagMaker.__accepted_reps).strip('[]')
            raise NotImplementedError(
                'Representation \'{}\' is unsupported. Accepted representations are {} .'.format(rep_str, accept_reps))
    
    def bty(self, dataset):
        '''
        Function to iterate over all molecules in a dataset
        and construct the bags that will be filled as well
        as what the size of each bag needs to be.
        Parameters
        ---------
        dataset: path
            path to all molecules in the dataset
        filetype: str
            molecule file type
        Returns
        -------
        bags: dict
            dict of all bags for the dataset
        bag_sizes: dict
            dict of size of the largest bags in the dataset
        '''
        # iterate through all of the molecules in the dataset
        #   and get the sizes of the largest bags
        bond_sizes = {}
        for molecule in glob.iglob("{}/*".format(dataset)):
            filetype = molecule.split('.')[-1]
            mol = next(pybel.readfile(filetype, molecule))

            # build bags
            bond_bag = {}

            # iterate through all bonds
            for bond in ob.OBMolBondIter(mol.OBMol):
                begin = atomType(mol, bond.GetBeginAtomIdx())
                end = atomType(mol, bond.GetEndAtomIdx())

                if (end < begin):
                    # swap them for lexographic order
                    begin, end = end, begin

                bond_type = encode[begin] + encode[end]
                if bond_type in bond_bag:
                    bond_bag[bond_type] += 1
                else:
                    bond_bag[bond_type] = 1

            # update bag_sizes with larger value
            bond_bag_key = list(bond_bag.keys())
            for i in range(len(bond_bag_key)):
                key = bond_bag_key[i]
                if key in bond_sizes:
                    if bond_bag[key] > bond_sizes[key]:
                        bond_sizes[key] = bond_bag[key]
                    else:
                        pass
                else:
                    bond_sizes[key] = bond_bag[key]

        # order bags alphabetically
        self.bag_sizes = OrderedDict(
            sorted(bond_sizes.items(), key=lambda t: t[0]))

        # make empty bags to fill
        self.bags = {}
        bag_keys = list(self.bag_sizes.keys())
        for i in range(len(bag_keys)):
            self.bags.update({bag_keys[i]: []})

    def baty(self, dataset):
        '''
        Function to iterate over all molecules in a dataset
        and construct the bags that will be filled as well
        as what the size of each bag needs to be.
        Parameters
        ---------
        dataset: path
            path to all molecules in the dataset
        filetype: str
            molecule file type
        Returns
        -------
        bags: dict
            dict of all bags for the dataset
        bag_sizes: dict
            dict of size of the largest bags in the dataset
        '''
        # iterate through all of the molecules in the dataset
        #   and get the sizes of the largest bags
        bond_sizes = {}
        angle_sizes = {}
        for molecule in glob.iglob("{}/*".format(dataset)):
            filetype = molecule.split('.')[-1]
            mol = next(pybel.readfile(filetype, molecule))

            # build bags
            bond_bag = {}
            angle_bag = {}

            # iterate through all bonds
            for bond in ob.OBMolBondIter(mol.OBMol):
                begin = atomType(mol, bond.GetBeginAtomIdx())
                end = atomType(mol, bond.GetEndAtomIdx())

                if (end < begin):
                    # swap them for lexographic order
                    begin, end = end, begin

                bond_type = encode[begin] + encode[end]
                if bond_type in bond_bag:
                    bond_bag[bond_type] += 1
                else:
                    bond_bag[bond_type] = 1

            # update bag_sizes with larger value
            bond_bag_key = list(bond_bag.keys())
            for i in range(len(bond_bag_key)):
                key = bond_bag_key[i]
                if key in bond_sizes:
                    if bond_bag[key] > bond_sizes[key]:
                        bond_sizes[key] = bond_bag[key]
                    else:
                        pass
                else:
                    bond_sizes[key] = bond_bag[key]

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

                if angle_type in angle_bag:
                    angle_bag[angle_type] += 1
                else:
                    angle_bag[angle_type] = 1

            # update bag_sizes with larger value
            angle_bag_key = list(angle_bag.keys())
            for i in range(len(angle_bag_key)):
                key = angle_bag_key[i]
                if key in angle_sizes:
                    if angle_bag[key] > angle_sizes[key]:
                        angle_sizes[key] = angle_bag[key]
                    else:
                        pass
                else:
                    angle_sizes[key] = angle_bag[key]

        self.bag_sizes = bond_sizes.copy()
        self.bag_sizes.update(angle_sizes)

        # order bags alphabetically
        self.bag_sizes = OrderedDict(
            sorted(self.bag_sizes.items(), key=lambda t: t[0]))


        # make empty bags to fill
        self.bags = {}
        bag_keys = list(self.bag_sizes.keys())
        for i in range(len(bag_keys)):
            self.bags.update({bag_keys[i]: []})

    def batty(self, dataset):
        '''
        Function to iterate over all molecules in a dataset
        and construct the bags that will be filled as well
        as what the size of each bag needs to be.
        Parameters
        ---------
        dataset: path
            path to all molecules in the dataset
        filetype: str
            molecule file type
        Returns
        -------
        bags: dict
            dict of all bags for the dataset
        bag_sizes: dict
            dict of size of the largest bags in the dataset
        '''
        # iterate through all of the molecules in the dataset
        #   and get the sizes of the largest bags
        bond_sizes = {}
        angle_sizes = {}
        torsion_sizes = {}
        for molecule in glob.iglob("{}/*".format(dataset)):
            filetype = molecule.split('.')[-1]
            mol = next(pybel.readfile(filetype, molecule))

            # build bags
            bond_bag = {}
            angle_bag = {}
            torsion_bag = {}

            # iterate through all bonds
            for bond in ob.OBMolBondIter(mol.OBMol):
                begin = atomType(mol, bond.GetBeginAtomIdx())
                end = atomType(mol, bond.GetEndAtomIdx())

                if (end < begin):
                    # swap them for lexographic order
                    begin, end = end, begin

                bond_type = encode[begin] + encode[end]
                if bond_type in bond_bag:
                    bond_bag[bond_type] += 1
                else:
                    bond_bag[bond_type] = 1

            # update bag_sizes with larger value
            bond_bag_key = list(bond_bag.keys())
            for i in range(len(bond_bag_key)):
                key = bond_bag_key[i]
                if key in bond_sizes:
                    if bond_bag[key] > bond_sizes[key]:
                        bond_sizes[key] = bond_bag[key]
                    else:
                        pass
                else:
                    bond_sizes[key] = bond_bag[key]

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

                if angle_type in angle_bag:
                    angle_bag[angle_type] += 1
                else:
                    angle_bag[angle_type] = 1

            # update bag_sizes with larger value
            angle_bag_key = list(angle_bag.keys())
            for i in range(len(angle_bag_key)):
                key = angle_bag_key[i]
                if key in angle_sizes:
                    if angle_bag[key] > angle_sizes[key]:
                        angle_sizes[key] = angle_bag[key]
                    else:
                        pass
                else:
                    angle_sizes[key] = angle_bag[key]

            # iterate through all torsions
            for torsion in ob.OBMolTorsionIter(mol.OBMol):
                a = (torsion[0] + 1)
                b = (torsion[1] + 1)
                c = (torsion[2] + 1)
                d = (torsion[3] + 1)

                aType = atomType(mol, a)
                bType = atomType(mol, b)
                cType = atomType(mol, c)
                dType = atomType(mol, d)

                # output in lexographic order
                if (aType < dType):
                    torsion_type = encode[aType] + encode[bType] \
                                    + encode[cType] + encode[dType]
                    if torsion_type in torsion_bag:
                        torsion_bag[torsion_type] += 1
                    else:
                        torsion_bag[torsion_type] = 1
                else:
                    torsion_type = encode[dType] + encode[cType] \
                                    + encode[bType] + encode[aType]
                    if torsion_type in torsion_bag:
                        torsion_bag[torsion_type] += 1
                    else:
                        torsion_bag[torsion_type] = 1

            # update bag_sizes with larger value
            torsion_bag_key = list(torsion_bag.keys())
            for i in range(len(torsion_bag_key)):
                key = torsion_bag_key[i]
                if key in torsion_sizes:
                    if torsion_bag[key] > torsion_sizes[key]:
                        torsion_sizes[key] = torsion_bag[key]
                    else:
                        pass
                else:
                    torsion_sizes[key] = torsion_bag[key]

        self.bag_sizes = bond_sizes.copy()
        self.bag_sizes.update(angle_sizes)
        self.bag_sizes.update(torsion_sizes)

        # order bags alphabetically
        self.bag_sizes = OrderedDict(
            sorted(self.bag_sizes.items(), key=lambda t: t[0]))


        # make empty bags to fill
        self.bags = {}
        bag_keys = list(self.bag_sizes.keys())
        for i in range(len(bag_keys)):
            self.bags.update({bag_keys[i]: []})

