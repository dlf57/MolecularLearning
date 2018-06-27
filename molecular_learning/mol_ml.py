#!/usr/bin/env python

import os
import pickle
import glob
import numpy as np
import pybel
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR
from sklearn.linear_model import BayesianRidge
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from itertools import repeat

# don't need to reimport openbabel
ob = pybel.ob

# Organization of all paths for the code
jobs_directory = '/home/dakota/Molecules/QM9/QM9-sdf/'  # molecule directory
jobs_file = '*.sdf'  # molecule files
jobs_file_format = 'sdf'  # molecule filetype

# import csv with properties
molec_csv = '/home/dakota/Molecules/QM9/qm9-rev2.csv'
csv_all = pd.read_csv(molec_csv)

# to get name split directory
start_split = jobs_directory.count('/')
end_split = start_split + 1


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


def molecular_info(mol, bond_len=200, angle_len=400, torsion_len=600, nb_len=1000, fill_val=-99999):
    # Pad lists to ensure same length
    # bond_len = 200
    # angle_len = 400
    # torsion_len = 600
    # nb_len = 1000
    # fill_val = -99999

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
    if (len(bonds) > bond_len):
        raise Exception('{} bonds. Change padding.'.format(len(bonds)))

    bonds.extend([fill_val] * (bond_len - len(bonds)))

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
    if (len(angles) > angle_len):
        raise Exception(
            '{} angles. Change padding.'.format(len(angles)))

    angles.extend([fill_val] * (angle_len - len(angles)))

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
    if (len(torsions) > torsion_len):
        raise Exception(
            '{} torsions. Change padding.'.format(len(torsions)))

    torsions.extend([fill_val] * (torsion_len - len(torsions)))

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
    # if (len(nb) > nb_len):
    #     raise Exception(
    #         '{} nonbonding. Change padding.'.format(len(nb)))
    if len(nb) < nb_len:
        nb.extend([fill_val] * (nb_len - len(nb)))
    else:
        nb = nb[:nb_len]

    # Create molecular descriptor
    n_rep = [*bonds, *angles, *torsions, *nb]

    return n_rep


row_list = []
# Read through all the files in the folder of this directory
for directory in sorted(glob.iglob(jobs_directory + jobs_file)):
    # name of the entry
    name = "/".join(directory.split('/')[start_split:end_split])
    short_name = name.split('.')[0]
    # name = name + '.xyz'

    try:
        mol = next(pybel.readfile(jobs_file_format, directory))

        # get the representation
        n_rep = molecular_info(mol)

        dict1 = {}
        dict1.update({'Descriptor': n_rep})
        dict1.update({'name': name})
        dict1.update({'Molecule': short_name})
        row_list.append(dict1)

    except StopIteration:
        pass

# dataframe of all of the bonds
df = pd.DataFrame(row_list, columns=['name', 'Molecule', 'Descriptor'])

# merge the two dataframes
combo = pd.merge(df, csv_all, on=['name'])

# molecular descriptors for ML
energy = np.asarray(list(combo[' dftE']), dtype=np.float)
molecule = np.asarray(list(combo['Descriptor']), dtype=np.float)

# test, train, split for data
X_train, X_test, y_train, y_test = train_test_split(
    molecule, energy, test_size=0.2)
# print(X_train.shape) # if array size errors print to see shape
# print(y_train.shape) # if array size errors print to see shape

# specify regressor and run
regr = RandomForestRegressor(n_jobs=-1)
regr.fit(X_train, y_train)
accuracy = regr.score(X_test, y_test)
r2 = ("%8.3f" % accuracy)
print('r^2 =', accuracy)
predicted = regr.predict(X_test)
mae = mean_absolute_error(y_test, predicted)
print('Mean Absolute Error =', mae)


# plot the actual versus the predicted
plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(predicted, y_test)
plt.annotate('$r^2$=' + str(r2), xy=(0.97, 0.10),
             xycoords='axes fraction', ha='right')
plt.annotate('$MAE$=' + str("%8.4f" % mae), xy=(0.97, 0.05),
             xycoords='axes fraction', ha='right')
plt.title('Comparison of DFT Energies', fontsize=16)
plt.xlabel('Predicted', fontsize=12)
plt.ylabel('Actual', fontsize=12)
plt.grid(False)
plt.savefig(save_figure)
plt.show()
plt.close
