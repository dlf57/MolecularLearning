'''
Function for returning molecular representations
'''

# Author: Dakota Folmsbee <dfolmsbee@gmail.com>
# License: GPLv2


import pybel

# don't need to reimport openbabel
ob = pybel.ob


def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()


# encoding dictionary
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


class Representation:
    def __init__(self, rep, n_bonds=50, n_angles=100, n_torsions=250, n_nb=400, fill_val=-99999):
        '''
        Representation class.

        Parameters
        ----------
        mol: pybel.Molecule object
            Molecule imported into pybel
        rep: str
            Defines which type of representation will be used
        n_bonds: int
            Sets size for bond list
        n_angles: int
            Sets size for angle list
        n_torsions: int
            Sets size for torsion list
        n_nb: int
            Sets size for nonbonding interaction list
        fill_val: int
            Arbitrary value to pad lists 
        '''
        rep_dict = {
            'BTY': self.BTy,
            'BATY': self.BATy,
            'BATTY': self.BATTy,
            'BATTYNB': self.BATTyNB
        }
        # self.mol = mol
        self.n_bonds = n_bonds
        self. n_angles = n_angles
        self.n_torsions = n_torsions
        self.n_nb = n_nb
        self.fill_val = fill_val
        self.f = rep_dict[rep]

    def BTy(self, mol):
        # iterate through all bonds
        bonds = []
        for bond in ob.OBMolBondIter(mol.OBMol):
            begin = atomType(mol, bond.GetBeginAtomIdx())
            end = atomType(mol, bond.GetEndAtomIdx())

            if (end < begin):
                # swap them for lexographic order
                begin, end = end, begin

            bond_type = encode[begin] + encode[end]
            bond_length = ("%8.4f" % (bond.GetLength()))
            bonds.append(float(bond_type))
            bonds.append(float(bond_length))

        # check to see if length does not exceed index
        if (len(bonds) > (self.n_bonds * 2)):
            raise Exception('{} bonds. Change padding.'.format(len(bonds)))

        bonds.extend([self.fill_val] * ((self.n_bonds * 2) - len(bonds)))

        # Create molecular descriptor
        b_rep = bonds

        return b_rep

    def BATy(self, mol):
        # iterate through all bonds
        bonds = Representation.BTy(self, mol)

        # iterate through all angles
        angles = []
        for angle in ob.OBMolAngleIter(mol.OBMol):
            a = (angle[0] + 1)
            b = mol.OBMol.GetAtom(angle[1] + 1)
            c = (angle[2] + 1)

            aType = atomType(mol, a)
            cType = atomType(mol, c)
            if (cType < aType):
                # swap them for lexographic order
                aType, cType = cType, aType

            angle_type = (encode[aType] + encode[b.GetType()]
                          + encode[cType])
            angle_angle = ("%8.3f" % (b.GetAngle(a, c)))
            angles.append(float(angle_type))
            angles.append(float(angle_angle))

        # check to see if length does not exceed index
        if (len(angles) > (self.n_angles * 2)):
            raise Exception(
                '{} angles. Change padding.'.format(len(angles)))

        angles.extend([self.fill_val] * ((self.n_angles * 2) - len(angles)))

        # Create molecular descriptor
        a_rep = [*bonds, *angles]

        return a_rep

    def BATTy(self, mol):
        # iterate through all bonds and angles
        bo_ang = Representation.BATy(self, mol)

        # iterate through all torsions
        torsions = []
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
                torsion_type = (encode[aType] + encode[bType]
                                + encode[cType] + encode[dType])
                torsion_angle = ("%8.3f" %
                                 ((mol.OBMol.GetTorsion(a, b, c, d))))
                torsions.append(float(torsion_type))
                torsions.append(float(torsion_angle))
            else:
                torsion_type = (encode[dType] + encode[cType]
                                + encode[bType] + encode[aType])
                torsion_angle = ("%8.3f" %
                                 ((mol.OBMol.GetTorsion(a, b, c, d))))
                torsions.append(float(torsion_type))
                torsions.append(float(torsion_angle))

        # check to see if length does not exceed index
        if (len(torsions) > (self.n_torsions * 2)):
            raise Exception(
                '{} torsions. Change padding.'.format(len(torsions)))

        torsions.extend([self.fill_val] *
                        ((self.n_torsions * 2) - len(torsions)))

        # Create molecular descriptor
        t_rep = [*bo_ang, *torsions]

        return t_rep

    def BATTyNB(self, mol):
        # iterate through all bonds, angles, and torsions
        bo_ang_tor = Representation.BATTy(self, mol)

        # iterate through all non-bonded
        nb = []
        for pair in ob.OBMolPairIter(mol.OBMol):
            (first, second) = pair
            begin_nb = atomType(mol, first)
            end_nb = atomType(mol, second)
            dist = mol.OBMol.GetAtom(first).GetDistance(second)
            if (end_nb[0] < begin_nb[0]):
                # swap them for lexographic order
                begin_nb, end_nb = end_nb, begin_nb

            nonbond_type = encode[begin_nb] + encode[end_nb]
            nonbond_distance = ("%8.4f" % (dist))
            nb.append(float(nonbond_type))
            nb.append(float(nonbond_distance))

        # check to see if length does not exceed index
        if (len(nb) > (self.n_nb * 2)):
            raise Exception(
                '{} nonbonding. Change padding.'.format(len(nb)))

        nb.extend([self.fill_val] * ((self.n_nb * 2) - len(nb)))

        # Create molecular descriptor
        nb_rep = [*bo_ang_tor, *nb]

        return nb_rep

    def rep(self, mol):
        return self.f(mol)
